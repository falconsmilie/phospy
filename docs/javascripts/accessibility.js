/* Accessible overflow behaviour for genuine wide data tables. */
(function () {
  "use strict";

  const observedWrappers = new WeakSet();

  function contextHeading(table) {
    const article = table.closest("article") || document;
    const headings = article.querySelectorAll("h1[id], h2[id], h3[id], h4[id]");
    let context = null;

    for (const heading of headings) {
      if (heading.compareDocumentPosition(table) & Node.DOCUMENT_POSITION_FOLLOWING) {
        context = heading;
      } else {
        break;
      }
    }
    return context;
  }

  function readableHeadingText(heading) {
    if (!heading) return "";
    const copy = heading.cloneNode(true);
    copy.querySelectorAll(".headerlink").forEach((link) => link.remove());
    return copy.textContent.trim();
  }

  function removeOverflowAccess(wrapper) {
    wrapper.classList.remove("phospy-table-scroll");
    wrapper.removeAttribute("role");
    wrapper.removeAttribute("tabindex");
    wrapper.removeAttribute("aria-label");
    wrapper.removeAttribute("aria-labelledby");
    wrapper.removeAttribute("aria-describedby");
    wrapper.previousElementSibling
      ?.matches(".phospy-table-scroll-hint") &&
      wrapper.previousElementSibling.remove();
  }

  function updateOverflowAccess(wrapper, table, index) {
    const overflows = wrapper.scrollWidth > wrapper.clientWidth + 1;
    if (!overflows) {
      removeOverflowAccess(wrapper);
      return;
    }

    const heading = contextHeading(table);
    const caption = table.caption?.textContent.trim();
    const contextText = caption || readableHeadingText(heading) || `Table ${index + 1}`;
    const hintId = `phospy-table-scroll-hint-${index + 1}`;

    wrapper.classList.add("phospy-table-scroll");
    wrapper.setAttribute("role", "region");
    wrapper.setAttribute("tabindex", "0");
    wrapper.setAttribute("aria-label", `Scrollable table: ${contextText}`);
    wrapper.removeAttribute("aria-labelledby");
    wrapper.setAttribute("aria-describedby", hintId);

    let hint = wrapper.previousElementSibling;
    if (!hint?.matches(".phospy-table-scroll-hint")) {
      hint = document.createElement("p");
      hint.className = "phospy-table-scroll-hint";
      wrapper.before(hint);
    }
    hint.id = hintId;
    hint.textContent = "Scrollable table. Scroll horizontally to see all columns.";
  }

  function enhanceTables() {
    document.querySelectorAll(".md-typeset__table").forEach((wrapper, index) => {
      const table = wrapper.querySelector(":scope > table") || wrapper.querySelector("table");
      if (!table) return;

      const columns = Math.max(
        0,
        ...Array.from(table.rows, (row) => row.cells.length),
      );
      wrapper.classList.toggle("phospy-table-wide", columns >= 4);
      wrapper.style.setProperty("--phospy-table-columns", String(columns));

      const update = () => updateOverflowAccess(wrapper, table, index);
      requestAnimationFrame(update);

      if (!observedWrappers.has(wrapper) && "ResizeObserver" in window) {
        new ResizeObserver(update).observe(wrapper);
        observedWrappers.add(wrapper);
      }
    });
  }

  if (typeof document$ !== "undefined") {
    document$.subscribe(enhanceTables);
  } else if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", enhanceTables, { once: true });
  } else {
    enhanceTables();
  }
})();
