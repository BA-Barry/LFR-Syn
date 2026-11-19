"""Shared inventory membership for CASP systems with different stock formats.

The shared stock artifact contains standard InChIKeys.  Syntheseus' default
inventory uses canonical SMILES, which is not strictly equivalent when multiple
SMILES collapse to one standard InChIKey.  This module provides the small
duck-typed inventory interface Syntheseus search algorithms require.  A separate
AiZynthFinder adapter is still required because its unsanitized runtime molecule
can lose an E/Z InChIKey layer even when the stock file itself is standardized.
"""

from __future__ import annotations

from dataclasses import dataclass
import re
from pathlib import Path
from typing import Dict, Protocol, Tuple

from rdkit import Chem


_INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")


class MoleculeLike(Protocol):
    smiles: str
    metadata: dict


class InchiKeyInventory:
    """Inventory compatible with Syntheseus' ``BaseMolInventory`` protocol.

    The class intentionally does not enumerate purchasable structures: an
    InChIKey is not reversible to SMILES.  Retrosynthesis search only needs
    ``is_purchasable`` and ``fill_metadata``.
    """

    def __init__(self, inchi_keys: set[str] | frozenset[str]):
        invalid = sorted(key for key in inchi_keys if not _INCHIKEY_RE.fullmatch(key))
        if invalid:
            preview = ", ".join(invalid[:3])
            raise ValueError(f"invalid standard InChIKey(s): {preview}")
        if not inchi_keys:
            raise ValueError("inventory must contain at least one InChIKey")
        self._inchi_keys = frozenset(inchi_keys)

    @classmethod
    def load_from_file(cls, path: str | Path) -> "InchiKeyInventory":
        lines = Path(path).read_text(encoding="utf-8").splitlines()
        if any(not line.strip() for line in lines):
            raise ValueError("inventory contains a blank InChIKey row")
        return cls({line.strip() for line in lines})

    def is_purchasable(self, mol: MoleculeLike) -> bool:
        return self.is_smiles_purchasable(mol.smiles)

    def is_smiles_purchasable(self, smiles: str) -> bool:
        """Query one explicit SMILES with the same standard-InChIKey rule."""

        if not isinstance(smiles, str):
            raise TypeError("inventory SMILES query must be a string")
        rd_mol = Chem.MolFromSmiles(smiles)
        if rd_mol is None:
            return False
        return Chem.MolToInchiKey(rd_mol) in self._inchi_keys

    def fill_metadata(self, mol: MoleculeLike) -> None:
        mol.metadata["is_purchasable"] = self.is_purchasable(mol)

    def __len__(self) -> int:
        return len(self._inchi_keys)

    @property
    def stock_inchikeys(self) -> frozenset[str]:
        return self._inchi_keys


@dataclass(frozen=True)
class StockRepresentation:
    source_row: int
    original_smiles: str
    canonical_smiles: str


@dataclass(frozen=True)
class EnumerableStockEntry:
    inchikey: str
    representations: Tuple[StockRepresentation, ...]


class EnumerableSmilesStock:
    """Enumerable view of the same logical standard-InChIKey inventory."""

    def __init__(self, entries: Tuple[EnumerableStockEntry, ...], raw_row_count: int):
        if not entries or raw_row_count < len(entries):
            raise ValueError("enumerable stock is empty or has inconsistent row counts")
        keys = [entry.inchikey for entry in entries]
        if len(set(keys)) != len(keys):
            raise ValueError("enumerable stock entries must have unique InChIKeys")
        self._entries = tuple(sorted(entries, key=lambda entry: entry.inchikey))
        self.raw_row_count = raw_row_count

    @classmethod
    def load_from_smiles_file(cls, path: str | Path) -> "EnumerableSmilesStock":
        lines = Path(path).read_text(encoding="utf-8").splitlines()
        if not lines:
            raise ValueError("SMILES stock must contain at least one row")
        if any(not line.strip() for line in lines):
            raise ValueError("SMILES stock contains a blank row")
        grouped: Dict[str, list[StockRepresentation]] = {}
        for index, original in enumerate(lines):
            smiles = original.strip()
            molecule = Chem.MolFromSmiles(smiles)
            if molecule is None:
                raise ValueError(f"invalid stock SMILES at row {index}: {smiles}")
            canonical = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
            inchikey = Chem.MolToInchiKey(molecule)
            grouped.setdefault(inchikey, []).append(
                StockRepresentation(index, smiles, canonical)
            )
        entries = tuple(
            EnumerableStockEntry(
                inchikey=inchikey,
                representations=tuple(
                    sorted(
                        representations,
                        key=lambda item: (item.canonical_smiles, item.source_row),
                    )
                ),
            )
            for inchikey, representations in grouped.items()
        )
        return cls(entries, len(lines))

    def assert_same_membership(self, inventory: InchiKeyInventory) -> None:
        observed = self.stock_inchikeys
        expected = inventory.stock_inchikeys
        if observed != expected:
            raise ValueError(
                "enumerable SMILES stock does not match logical InChIKey inventory: "
                f"missing={len(expected - observed)},extra={len(observed - expected)}"
            )

    @property
    def entries(self) -> Tuple[EnumerableStockEntry, ...]:
        return self._entries

    @property
    def stock_inchikeys(self) -> frozenset[str]:
        return frozenset(entry.inchikey for entry in self._entries)

    def __len__(self) -> int:
        return len(self._entries)
