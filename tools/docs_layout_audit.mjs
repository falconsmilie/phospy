#!/usr/bin/env node

import { spawn } from "node:child_process";
import { createServer } from "node:http";
import { mkdtemp, readdir, readFile, rm, writeFile } from "node:fs/promises";
import { tmpdir } from "node:os";
import { extname, join, relative, resolve, sep } from "node:path";

const MIME_TYPES = {
  ".css": "text/css; charset=utf-8",
  ".html": "text/html; charset=utf-8",
  ".js": "text/javascript; charset=utf-8",
  ".json": "application/json; charset=utf-8",
  ".png": "image/png",
  ".svg": "image/svg+xml",
  ".webmanifest": "application/manifest+json",
};

const siteRoot = resolve(process.argv[2] ?? "site");
const reportPath = resolve(process.argv[3] ?? "build/docs-layout-audit.json");
const screenshotRoot = process.argv[4] ? resolve(process.argv[4]) : null;
const chromePath =
  process.env.CHROME_PATH ??
  "C:\\Program Files\\Google\\Chrome\\Application\\chrome.exe";

async function htmlRoutes(root) {
  const routes = [];
  async function visit(directory) {
    for (const entry of await readdir(directory, { withFileTypes: true })) {
      const path = join(directory, entry.name);
      if (entry.isDirectory()) await visit(path);
      if (!entry.isFile() || entry.name !== "index.html") continue;
      const route = relative(root, path)
        .split(sep)
        .slice(0, -1)
        .join("/");
      routes.push(route ? `/${route}/` : "/");
    }
  }
  await visit(root);
  return routes.sort();
}

async function startStaticServer(root) {
  const server = createServer(async (request, response) => {
    try {
      const requestUrl = new URL(request.url, "http://127.0.0.1");
      let pathname = decodeURIComponent(requestUrl.pathname);
      if (pathname.endsWith("/")) pathname += "index.html";
      const target = resolve(root, `.${pathname}`);
      if (!target.startsWith(`${root}${sep}`) && target !== root) {
        response.writeHead(403).end("Forbidden");
        return;
      }
      const body = await readFile(target);
      response.writeHead(200, {
        "content-type": MIME_TYPES[extname(target)] ?? "application/octet-stream",
      });
      response.end(body);
    } catch {
      response.writeHead(404).end("Not found");
    }
  });
  await new Promise((resolveReady) => server.listen(0, "127.0.0.1", resolveReady));
  const { port } = server.address();
  return { server, origin: `http://127.0.0.1:${port}` };
}

async function startChrome() {
  const profile = await mkdtemp(join(tmpdir(), "phospy-docs-chrome-"));
  const child = spawn(
    chromePath,
    [
      "--headless=new",
      "--disable-gpu",
      "--disable-software-rasterizer",
      "--disable-extensions",
      "--hide-scrollbars=false",
      "--no-first-run",
      "--no-default-browser-check",
      "--no-sandbox",
      "--remote-debugging-port=0",
      "--remote-allow-origins=*",
      `--user-data-dir=${profile}`,
      "about:blank",
    ],
    { stdio: ["ignore", "ignore", "pipe"], windowsHide: true },
  );
  const browserWebSocket = await new Promise((resolveReady, reject) => {
    let stderr = "";
    const timer = setTimeout(() => reject(new Error("Chrome did not start")), 15000);
    child.stderr.on("data", (chunk) => {
      stderr += chunk.toString();
      const match = stderr.match(/DevTools listening on (ws:\/\/[^\s]+)/);
      if (!match) return;
      clearTimeout(timer);
      resolveReady(match[1]);
    });
    child.once("exit", (code) => {
      clearTimeout(timer);
      reject(new Error(`Chrome exited before DevTools was ready (${code})`));
    });
  });
  const port = new URL(browserWebSocket).port;
  const pages = await fetch(`http://127.0.0.1:${port}/json/list`).then((response) =>
    response.json(),
  );
  const page = pages.find((target) => target.type === "page");
  if (!page) throw new Error("Chrome did not expose a page target");
  return {
    child,
    profile,
    pageWebSocket: page.webSocketDebuggerUrl,
    async close() {
      child.kill();
      if (child.exitCode === null) {
        await Promise.race([
          new Promise((resolveExit) => child.once("exit", resolveExit)),
          new Promise((resolveTimeout) => setTimeout(resolveTimeout, 3000)),
        ]);
      }
      try {
        await rm(profile, { recursive: true, force: true, maxRetries: 5, retryDelay: 100 });
      } catch (error) {
        console.warn(`Could not remove temporary Chrome profile: ${error.message}`);
      }
    },
  };
}

class DevTools {
  constructor(url) {
    this.nextId = 1;
    this.pending = new Map();
    this.events = new Map();
    this.socket = new WebSocket(url);
    this.socket.addEventListener("message", ({ data }) => {
      const message = JSON.parse(data);
      if (message.id) {
        const pending = this.pending.get(message.id);
        this.pending.delete(message.id);
        if (message.error) pending.reject(new Error(message.error.message));
        else pending.resolve(message.result);
        return;
      }
      const listeners = this.events.get(message.method) ?? [];
      for (const listener of listeners.splice(0)) listener(message.params);
    });
    this.socket.addEventListener("close", () => {
      for (const pending of this.pending.values()) {
        pending.reject(new Error("Chrome DevTools connection closed"));
      }
      this.pending.clear();
    });
  }

  async open() {
    if (this.socket.readyState === WebSocket.OPEN) return;
    await new Promise((resolveReady, reject) => {
      this.socket.addEventListener("open", resolveReady, { once: true });
      this.socket.addEventListener("error", reject, { once: true });
    });
  }

  send(method, params = {}) {
    const id = this.nextId++;
    return new Promise((resolveResult, reject) => {
      this.pending.set(id, { resolve: resolveResult, reject });
      this.socket.send(JSON.stringify({ id, method, params }));
    });
  }

  once(method) {
    return new Promise((resolveEvent) => {
      const listeners = this.events.get(method) ?? [];
      listeners.push(resolveEvent);
      this.events.set(method, listeners);
    });
  }

  close() {
    this.socket.close();
  }
}

const auditExpression = `(() => {
  const root = document.documentElement;
  const width = root.clientWidth;
  const visible = (element) => {
    const style = getComputedStyle(element);
    const rect = element.getBoundingClientRect();
    return style.display !== "none" && style.visibility !== "hidden" && style.opacity !== "0" && rect.width > 0 && rect.height > 0;
  };
  const selector = (element) => {
    if (element.id) return "#" + CSS.escape(element.id);
    const classes = [...element.classList].slice(0, 3).map((name) => "." + CSS.escape(name)).join("");
    return element.tagName.toLowerCase() + classes;
  };
  const offenders = [...document.body.querySelectorAll("*")]
    .filter(visible)
    .map((element) => ({ element, rect: element.getBoundingClientRect() }))
    .filter(({ element }) => getComputedStyle(element).pointerEvents !== "none" && !element.disabled)
    .filter(({ rect }) => rect.right > width + 1 || rect.left < -1)
    .filter(({ element }) => !element.closest(".md-typeset__table, pre, .highlight"))
    .slice(0, 12)
    .map(({ element, rect }) => ({
      selector: selector(element),
      left: Math.round(rect.left),
      right: Math.round(rect.right),
      width: Math.round(rect.width),
      whiteSpace: getComputedStyle(element).whiteSpace,
      minWidth: getComputedStyle(element).minWidth,
      overflowX: getComputedStyle(element).overflowX,
      text: (element.textContent || "").trim().replace(/\s+/g, " ").slice(0, 100),
    }));
  const tables = [...document.querySelectorAll("table")].map((table, index) => {
    const wrapper = table.closest(".md-typeset__table");
    const cells = [...table.querySelectorAll("th, td")];
    const widestCell = cells
      .map((cell) => ({
        text: (cell.textContent || "").trim().replace(/\s+/g, " ").slice(0, 120),
        width: Math.round(cell.getBoundingClientRect().width),
        whiteSpace: getComputedStyle(cell).whiteSpace,
        overflowWrap: getComputedStyle(cell).overflowWrap,
      }))
      .sort((left, right) => right.width - left.width)[0];
    return {
      index,
      columns: table.tHead?.rows[0]?.cells.length ?? table.rows[0]?.cells.length ?? 0,
      rows: table.rows.length,
      tableWidth: Math.round(table.getBoundingClientRect().width),
      wrapperWidth: wrapper ? Math.round(wrapper.getBoundingClientRect().width) : null,
      scrollWidth: wrapper?.scrollWidth ?? null,
      clientWidth: wrapper?.clientWidth ?? null,
      overflowX: wrapper ? getComputedStyle(wrapper).overflowX : null,
      focusableWrapper: wrapper?.tabIndex === 0,
      wrapperName: wrapper?.getAttribute("aria-label") ?? wrapper?.getAttribute("aria-labelledby"),
      caption: table.caption?.textContent.trim() ?? null,
      widestCell,
    };
  });
  const interactiveWithoutName = [...document.querySelectorAll("a, button, input, select, textarea")]
    .filter(visible)
    .filter((element) => {
      const labelled = element.getAttribute("aria-label") || element.getAttribute("aria-labelledby") || element.getAttribute("title");
      const text = (element.textContent || "").trim();
      const value = element.value || element.getAttribute("placeholder") || element.getAttribute("alt");
      return !(labelled || text || value);
    })
    .map(selector);
  const smallPointerTargets = [...document.querySelectorAll("button, input, select, textarea, summary, label[for], .md-header__button")]
    .filter(visible)
    .map((element) => ({ element, rect: element.getBoundingClientRect() }))
    .filter(({ rect }) => rect.width < 24 || rect.height < 24)
    .map(({ element, rect }) => ({
      selector: selector(element),
      width: Math.round(rect.width),
      height: Math.round(rect.height),
    }));
  const headings = [...document.querySelectorAll("h1, h2, h3, h4, h5, h6")].map((heading) => ({
    level: Number(heading.tagName.slice(1)),
    text: heading.textContent.trim(),
  }));
  const skippedHeadingLevel = headings.some((heading, index) => index > 0 && heading.level > headings[index - 1].level + 1);
  const skipLink = document.querySelector('a[href="#main-content"], a[href^="#"]');
  return {
    title: document.title,
    viewportWidth: width,
    scrollWidth: root.scrollWidth,
    pageOverflow: root.scrollWidth > width + 1,
    offenders,
    tables,
    landmarks: {
      header: document.querySelectorAll("header").length,
      nav: document.querySelectorAll("nav").length,
      main: document.querySelectorAll("main").length,
      footer: document.querySelectorAll("footer").length,
    },
    h1Count: document.querySelectorAll("h1").length,
    skippedHeadingLevel,
    skipLink: skipLink ? { text: skipLink.textContent.trim(), href: skipLink.getAttribute("href") } : null,
    interactiveWithoutName,
    smallPointerTargets,
  };
})()`;

async function auditPage(devtools, url, width) {
  await devtools.send("Emulation.setDeviceMetricsOverride", {
    width,
    height: 900,
    deviceScaleFactor: 1,
    mobile: false,
  });
  await devtools.send("Page.navigate", { url });
  let ready = false;
  for (let attempt = 0; attempt < 1200; attempt += 1) {
    try {
      const state = await devtools.send("Runtime.evaluate", {
        expression:
          "({ href: location.href, ready: document.readyState, hasBody: document.body !== null })",
        returnByValue: true,
      });
      if (
        state.result.value?.href === url &&
        state.result.value?.hasBody &&
        state.result.value?.ready !== "loading"
      ) {
        ready = true;
        break;
      }
    } catch {
      // A navigation can replace the execution context between request and response.
    }
    await new Promise((resolveReady) => setTimeout(resolveReady, 30));
  }
  if (!ready) throw new Error(`Timed out loading ${url}`);
  await new Promise((resolveReady) => setTimeout(resolveReady, 80));
  const result = await devtools.send("Runtime.evaluate", {
    expression: auditExpression,
    returnByValue: true,
  });
  if (!result.result.value) {
    throw new Error(
      `Could not audit ${url} at ${width}px: ${JSON.stringify(result.exceptionDetails ?? result.result)}`,
    );
  }
  return result.result.value;
}

async function auditVariant(devtools, origin, variant) {
  await auditPage(devtools, `${origin}${variant.route}`, variant.width);
  await devtools.send("Runtime.evaluate", {
    expression: `(() => {
      const options = ${JSON.stringify({
        scheme: variant.scheme,
        textScale: variant.textScale,
        textSpacing: variant.textSpacing,
      })};
      if (options.scheme) document.body.setAttribute("data-md-color-scheme", options.scheme);
      if (options.textScale) document.documentElement.style.fontSize = options.textScale;
      if (options.textSpacing) {
        const style = document.createElement("style");
        style.id = "phospy-text-spacing-audit";
        style.textContent = "* { line-height: 1.5 !important; letter-spacing: 0.12em !important; word-spacing: 0.16em !important; } p { margin-bottom: 2em !important; }";
        document.head.append(style);
      }
    })()`,
  });
  await new Promise((resolveReady) => setTimeout(resolveReady, 120));
  const result = await devtools.send("Runtime.evaluate", {
    expression: auditExpression,
    returnByValue: true,
  });
  return { ...variant, ...result.result.value };
}

async function keyboardChecks(devtools, origin) {
  await auditPage(devtools, `${origin}/api/differential-analysis/`, 320);
  const focusedTable = await devtools.send("Runtime.evaluate", {
    expression: `(() => {
      const wrapper = document.querySelector(".phospy-table-scroll");
      if (!wrapper) return null;
      wrapper.scrollLeft = 0;
      wrapper.focus();
      return {
        active: document.activeElement === wrapper,
        role: wrapper.getAttribute("role"),
        tabindex: wrapper.getAttribute("tabindex"),
        name: wrapper.getAttribute("aria-label") || wrapper.getAttribute("aria-labelledby"),
        description: wrapper.getAttribute("aria-describedby"),
        outlineWidth: getComputedStyle(wrapper).outlineWidth,
        before: wrapper.scrollLeft,
      };
    })()`,
    returnByValue: true,
  });
  await devtools.send("Input.dispatchKeyEvent", {
    type: "keyDown",
    key: "ArrowRight",
    code: "ArrowRight",
    windowsVirtualKeyCode: 39,
    nativeVirtualKeyCode: 39,
  });
  await devtools.send("Input.dispatchKeyEvent", {
    type: "keyUp",
    key: "ArrowRight",
    code: "ArrowRight",
    windowsVirtualKeyCode: 39,
    nativeVirtualKeyCode: 39,
  });
  await new Promise((resolveReady) => setTimeout(resolveReady, 150));
  const scrolledTable = await devtools.send("Runtime.evaluate", {
    expression: `(() => {
      const wrapper = document.querySelector(".phospy-table-scroll");
      return wrapper ? { after: wrapper.scrollLeft, max: wrapper.scrollWidth - wrapper.clientWidth } : null;
    })()`,
    returnByValue: true,
  });

  await devtools.send("Accessibility.enable");
  const accessibilityTree = await devtools.send("Accessibility.getFullAXTree");
  const roleCounts = {};
  const regions = [];
  for (const node of accessibilityTree.nodes) {
    const role = node.role?.value;
    if (["table", "row", "columnheader", "cell", "region"].includes(role)) {
      roleCounts[role] = (roleCounts[role] ?? 0) + 1;
    }
    if (role === "region") regions.push(node.name?.value ?? "");
  }

  await auditPage(devtools, `${origin}/`, 1280);
  await devtools.send("Runtime.evaluate", {
    expression: "document.body.focus(); document.activeElement?.blur();",
  });
  const tabStops = [];
  for (let index = 0; index < 7; index += 1) {
    await devtools.send("Input.dispatchKeyEvent", {
      type: "keyDown",
      key: "Tab",
      code: "Tab",
      windowsVirtualKeyCode: 9,
      nativeVirtualKeyCode: 9,
    });
    await devtools.send("Input.dispatchKeyEvent", {
      type: "keyUp",
      key: "Tab",
      code: "Tab",
      windowsVirtualKeyCode: 9,
      nativeVirtualKeyCode: 9,
    });
    const active = await devtools.send("Runtime.evaluate", {
      expression: `(() => {
        const element = document.activeElement;
        const rect = element.getBoundingClientRect();
        return {
          tag: element.tagName.toLowerCase(),
          className: element.className,
          label: element.getAttribute("aria-label") || element.getAttribute("title") || element.textContent.trim().replace(/\\s+/g, " ").slice(0, 80),
          visibleInViewport: rect.bottom > 0 && rect.top < innerHeight && rect.right > 0 && rect.left < innerWidth,
          outlineWidth: getComputedStyle(element).outlineWidth,
        };
      })()`,
      returnByValue: true,
    });
    tabStops.push(active.result.value);
  }

  const zoomMetrics = async () => {
    const result = await devtools.send("Runtime.evaluate", {
      expression:
        "({ innerWidth, outerWidth, devicePixelRatio, visualScale: visualViewport?.scale, scrollWidth: document.documentElement.scrollWidth })",
      returnByValue: true,
    });
    return result.result.value;
  };
  const browserZoomBefore = await zoomMetrics();
  for (let index = 0; index < 8; index += 1) {
    await devtools.send("Input.dispatchKeyEvent", {
      type: "keyDown",
      key: "+",
      code: "Equal",
      modifiers: 2,
      windowsVirtualKeyCode: 187,
      nativeVirtualKeyCode: 187,
    });
    await devtools.send("Input.dispatchKeyEvent", {
      type: "keyUp",
      key: "+",
      code: "Equal",
      modifiers: 2,
      windowsVirtualKeyCode: 187,
      nativeVirtualKeyCode: 187,
    });
  }
  await new Promise((resolveReady) => setTimeout(resolveReady, 250));
  const browserZoomAfter = await zoomMetrics();
  await devtools.send("Input.dispatchKeyEvent", {
    type: "keyDown",
    key: "0",
    code: "Digit0",
    modifiers: 2,
    windowsVirtualKeyCode: 48,
    nativeVirtualKeyCode: 48,
  });
  await devtools.send("Input.dispatchKeyEvent", {
    type: "keyUp",
    key: "0",
    code: "Digit0",
    modifiers: 2,
    windowsVirtualKeyCode: 48,
    nativeVirtualKeyCode: 48,
  });

  return {
    tableRegion: { ...focusedTable.result.value, ...scrolledTable.result.value },
    accessibilityTree: { roleCounts, regions },
    sharedControlTabStops: tabStops,
    browserZoomAttempt: { before: browserZoomBefore, after: browserZoomAfter },
  };
}

async function main() {
  const routes = await htmlRoutes(siteRoot);
  console.log(`Found ${routes.length} rendered routes.`);
  const { server, origin } = await startStaticServer(siteRoot);
  console.log(`Serving the rendered site at ${origin}.`);
  const chrome = await startChrome();
  console.log("Chrome DevTools connection is ready.");
  const devtools = new DevTools(chrome.pageWebSocket);
  await devtools.open();
  console.log("Enabling page inspection.");
  await devtools.send("Page.enable");
  await devtools.send("Runtime.enable");
  console.log("Page inspection is enabled.");
  const results = [];
  const variants = [];
  let interactions = null;
  try {
    for (const route of routes) {
      for (const width of [1280, 320]) {
        results.push({ route, ...(await auditPage(devtools, `${origin}${route}`, width)) });
      }
      if (results.length % 20 === 0) console.log(`Measured ${results.length} layouts...`);
    }
    for (const variant of [
      { name: "light-home", route: "/", width: 1280, scheme: "default" },
      { name: "light-guide-mobile", route: "/api/guide/", width: 320, scheme: "default" },
      { name: "light-reference-mobile", route: "/scientific-coverage/", width: 320, scheme: "default" },
      { name: "text-200-home", route: "/", width: 1280, textScale: "200%" },
      { name: "text-200-tutorial", route: "/quickstart/", width: 1280, textScale: "200%" },
      { name: "text-200-workflow", route: "/api/differential-analysis/", width: 1280, textScale: "200%" },
      { name: "text-200-reference", route: "/scientific-coverage/", width: 1280, textScale: "200%" },
      { name: "text-spacing-home", route: "/", width: 320, textSpacing: true },
      { name: "text-spacing-tutorial", route: "/quickstart/", width: 320, textSpacing: true },
      { name: "text-spacing-workflow", route: "/api/differential-analysis/", width: 320, textSpacing: true },
      { name: "text-spacing-reference", route: "/scientific-coverage/", width: 320, textSpacing: true },
    ]) {
      variants.push(await auditVariant(devtools, origin, variant));
    }
    interactions = await keyboardChecks(devtools, origin);
    if (screenshotRoot) {
      await import("node:fs/promises").then(({ mkdir }) => mkdir(screenshotRoot, { recursive: true }));
      for (const [name, route, width, height, scheme] of [
        ["home-desktop", "/", 1280, 900, "slate"],
        ["home-light-desktop", "/", 1280, 900, "default"],
        ["quickstart-mobile", "/quickstart/", 320, 800, "slate"],
        ["differential-tablet", "/api/differential-analysis/", 768, 900, "slate"],
        ["api-mobile", "/api/guide/", 320, 800, "slate"],
      ]) {
        await devtools.send("Emulation.setDeviceMetricsOverride", {
          width,
          height,
          deviceScaleFactor: 1,
          mobile: false,
        });
        await devtools.send("Page.navigate", { url: `${origin}${route}` });
        for (let attempt = 0; attempt < 60; attempt += 1) {
          const state = await devtools.send("Runtime.evaluate", {
            expression: "document.readyState",
            returnByValue: true,
          });
          if (state.result.value === "complete") break;
          await new Promise((resolveReady) => setTimeout(resolveReady, 50));
        }
        await devtools.send("Runtime.evaluate", {
          expression: `document.body.setAttribute("data-md-color-scheme", ${JSON.stringify(scheme)})`,
        });
        await new Promise((resolveReady) => setTimeout(resolveReady, 80));
        const capture = await devtools.send("Page.captureScreenshot", {
          format: "png",
          captureBeyondViewport: false,
        });
        await writeFile(join(screenshotRoot, `${name}.png`), capture.data, "base64");
      }
    }
  } finally {
    devtools.close();
    await chrome.close();
    await new Promise((resolveClosed) => server.close(resolveClosed));
  }
  const report = {
    generatedAt: new Date().toISOString(),
    siteRoot,
    routes: routes.length,
    results,
    variants,
    interactions,
  };
  await writeFile(reportPath, `${JSON.stringify(report, null, 2)}\n`, "utf8");
  const overflow = results.filter((result) => result.pageOverflow);
  const unnamed = results.filter((result) => result.interactiveWithoutName?.length);
  console.log(`Audited ${routes.length} routes at 1280px and 320px.`);
  console.log(`Page-wide overflow results: ${overflow.length}`);
  for (const result of overflow.slice(0, 20)) {
    console.log(`  ${result.viewportWidth}px ${result.route}: ${result.scrollWidth}px`);
  }
  console.log(`Results with unnamed interactive controls: ${unnamed.length}`);
  console.log(`Variant overflow results: ${variants.filter((result) => result.pageOverflow).length}`);
  console.log(`Keyboard-scrolled table position: ${interactions?.tableRegion?.after ?? "unavailable"}`);
  console.log(`Report: ${reportPath}`);
}

await main();
