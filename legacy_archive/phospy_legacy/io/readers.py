from __future__ import annotations

import json
import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    from ..api.contracts import (
        SimpleKinaseWorkflowBundleMetadata,
        WorkflowOutputInventoryItem,
    )

from ..internal.constants import (
    SIMPLE_KINASE_WORKFLOW_BUNDLE_FORMAT,
    SIMPLE_KINASE_WORKFLOW_RESULT_TYPE,
    WORKFLOW_OUTPUT_BUNDLE_MANIFEST_FILENAME,
)
from ..validation.schema.tables import (
    PhosphoInputSchema,
    PredMatSchema,
    TotalInputSchema,
)

DEFAULT_TEXT_ENCODING = "utf-8"


def clean_columns(columns: Iterable[str]) -> list[str]:
    cleaned: list[str] = []
    for col in columns:
        value = col.strip().lower()
        value = re.sub(r"[^0-9a-zA-Z]+", "_", value)
        value = re.sub(r"_+", "_", value)
        value = value.strip("_")
        cleaned.append(value)
    return cleaned


def default_text_encoding(path: str | Path | None = None) -> str:
    """Return the package default text encoding.

    The loader does not infer encodings from file contents. Callers should pass
    an explicit encoding when they need something other than the package
    default. The optional ``path`` argument is accepted for API convenience.
    """

    _ = path
    return DEFAULT_TEXT_ENCODING


def read_table_raw(
    path: str | Path,
    *,
    sep: str = "\t",
    encoding: str | None = None,
    index_col: int | str | None = None,
    usecols: Sequence[str | int] | None = None,
) -> pd.DataFrame:
    resolved_encoding = encoding or DEFAULT_TEXT_ENCODING
    return pd.read_csv(
        path,
        sep=sep,
        encoding=resolved_encoding,
        low_memory=False,
        index_col=index_col,
        usecols=usecols,
    )


def read_table(
    path: str | Path,
    encoding: str | None = None,
    *,
    usecols: Sequence[str | int] | None = None,
) -> pd.DataFrame:
    frame = read_table_raw(path, encoding=encoding, usecols=usecols)
    return clean_table_columns(frame)


def load_total_table(
    path: str | Path,
    *,
    encoding: str | None = None,
    usecols: Sequence[str | int] | None = None,
) -> pd.DataFrame:
    frame = clean_table_columns(
        read_table_raw(path, encoding=encoding, usecols=usecols)
    )
    return TotalInputSchema.validate(frame, context=f"total input table ({path})")


def load_phospho_table(
    path: str | Path,
    *,
    encoding: str | None = None,
    usecols: Sequence[str | int] | None = None,
) -> pd.DataFrame:
    frame = clean_table_columns(
        read_table_raw(path, encoding=encoding, usecols=usecols)
    )
    return PhosphoInputSchema.validate(frame, context=f"phospho input table ({path})")


def load_pred_mat(
    path: str | Path,
    *,
    encoding: str | None = None,
    usecols: Sequence[str | int] | None = None,
) -> pd.DataFrame:
    frame = read_table_raw(
        path,
        sep=",",
        encoding=encoding,
        index_col=0,
        usecols=usecols,
    )
    frame.index = frame.index.map(str)
    frame.columns = [str(column).strip() for column in frame.columns]
    return PredMatSchema.validate(frame, context=f"pred_mat ({path})")


def clean_table_columns(frame: pd.DataFrame) -> pd.DataFrame:
    frame.columns = clean_columns(str(column) for column in frame.columns)
    return frame


@dataclass(slots=True)
class SimpleKinaseWorkflowOutputBundle:
    """Loaded workflow output bundle payload."""

    metadata: SimpleKinaseWorkflowBundleMetadata
    tables: dict[str, pd.DataFrame | pd.Series]

    def get_table(self, table_id: str) -> pd.DataFrame | pd.Series:
        return self.tables[table_id]


@dataclass(frozen=True, slots=True)
class SimpleKinaseWorkflowBundleReader:
    """Read saved ``SimpleKinaseWorkflowResult`` output bundles."""

    def read_metadata(
        self,
        bundle_dir: str | Path,
    ) -> SimpleKinaseWorkflowBundleMetadata:
        from ..api.contracts import SimpleKinaseWorkflowBundleMetadata

        root_dir = Path(bundle_dir)
        manifest_path = root_dir / WORKFLOW_OUTPUT_BUNDLE_MANIFEST_FILENAME
        raw_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if not isinstance(raw_manifest, Mapping):
            msg = f"Bundle manifest {manifest_path} must contain a JSON object."
            raise ValueError(msg)
        metadata = SimpleKinaseWorkflowBundleMetadata.from_mapping(raw_manifest)
        self._validate_metadata(metadata, manifest_path=manifest_path)
        return metadata

    def read(
        self,
        bundle_dir: str | Path,
        *,
        table_ids: Sequence[str] | None = None,
    ) -> SimpleKinaseWorkflowOutputBundle:
        root_dir = Path(bundle_dir)
        metadata = self.read_metadata(root_dir)
        inventory = {item.table_id: item for item in metadata.output_inventory}

        selected_table_ids = (
            tuple(inventory)
            if table_ids is None
            else tuple(str(item) for item in table_ids)
        )
        unknown_table_ids = sorted(
            table_id for table_id in selected_table_ids if table_id not in inventory
        )
        if unknown_table_ids:
            available = ", ".join(sorted(inventory))
            requested = ", ".join(unknown_table_ids)
            msg = (
                f"Unknown output bundle table IDs: {requested}. "
                f"Available table IDs: {available}."
            )
            raise ValueError(msg)

        tables: dict[str, pd.DataFrame | pd.Series] = {}
        for table_id in selected_table_ids:
            item = inventory[table_id]
            table_path = root_dir / item.path
            tables[table_id] = self._read_inventory_item(item, table_path=table_path)

        return SimpleKinaseWorkflowOutputBundle(
            metadata=metadata,
            tables=tables,
        )

    @staticmethod
    def _read_inventory_item(
        item: WorkflowOutputInventoryItem,
        *,
        table_path: Path,
    ) -> pd.DataFrame | pd.Series:
        frame = pd.read_csv(table_path, index_col=0)
        if isinstance(frame.index.name, str) and frame.index.name.startswith("Unnamed"):
            frame.index.name = None

        if item.value_type == "dataframe":
            return frame

        if frame.shape[1] != 1:
            msg = f"Series table {item.table_id!r} in {table_path} must contain exactly one column."
            raise ValueError(msg)
        series = frame.iloc[:, 0]
        if isinstance(series.name, str) and series.name.startswith("Unnamed"):
            series.name = None
        return series

    @staticmethod
    def _validate_metadata(
        metadata: SimpleKinaseWorkflowBundleMetadata,
        *,
        manifest_path: Path,
    ) -> None:
        if metadata.workflow_type != SIMPLE_KINASE_WORKFLOW_RESULT_TYPE:
            msg = (
                f"Bundle manifest {manifest_path} has unsupported workflow type "
                f"{metadata.workflow_type!r}."
            )
            raise ValueError(msg)
        if metadata.bundle_format != SIMPLE_KINASE_WORKFLOW_BUNDLE_FORMAT:
            msg = (
                f"Bundle manifest {manifest_path} has unsupported bundle format "
                f"{metadata.bundle_format!r}."
            )
            raise ValueError(msg)


def load_simple_kinase_workflow_output_bundle_metadata(
    bundle_dir: str | Path,
) -> SimpleKinaseWorkflowBundleMetadata:
    return SimpleKinaseWorkflowBundleReader().read_metadata(bundle_dir)


def load_simple_kinase_workflow_output_bundle(
    bundle_dir: str | Path,
    *,
    table_ids: Sequence[str] | None = None,
) -> SimpleKinaseWorkflowOutputBundle:
    return SimpleKinaseWorkflowBundleReader().read(
        bundle_dir,
        table_ids=table_ids,
    )


__all__ = [
    "DEFAULT_TEXT_ENCODING",
    "SimpleKinaseWorkflowBundleReader",
    "SimpleKinaseWorkflowOutputBundle",
    "clean_columns",
    "clean_table_columns",
    "default_text_encoding",
    "load_simple_kinase_workflow_output_bundle",
    "load_simple_kinase_workflow_output_bundle_metadata",
    "load_phospho_table",
    "load_pred_mat",
    "load_total_table",
    "read_table",
    "read_table_raw",
]
