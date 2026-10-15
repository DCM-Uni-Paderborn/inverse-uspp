#!/usr/bin/env python3
"""UPF numeric comparison and byte-preserving XML helpers.

Derived from CP2K's add_upf_full_waves.py. Reference supplementation and
route selection are implemented in prepare_upf.py.
"""
# SPDX-License-Identifier: GPL-2.0-or-later

import math
import sys
from xml.parsers import expat


def real(text):
    value = float(text.replace("D", "E").replace("d", "e"))
    if not math.isfinite(value):
        raise ValueError("Nonfinite UPF numeric value")
    return value


def values(node):
    if list(node):
        return []
    result = [real(word) for word in (node.text or "").split()]
    if "size" in node.attrib and int(node.attrib["size"]) != len(result):
        raise ValueError(f"Incorrect array size in {node.tag}")
    return result


def flag(header, name):
    value = header.get(name, "false").strip().lower().strip(".")
    if value not in ("true", "false", "t", "f", "1", "0"):
        raise ValueError(f"Invalid logical attribute {name}")
    return value in ("true", "t", "1")


def normalized(text):
    try:
        value = float(text.replace("D", "E").replace("d", "e"))
    except ValueError:
        return " ".join(text.split()).lower()
    if not math.isfinite(value):
        raise ValueError("Nonfinite UPF numeric attribute")
    return value


def children(node):
    result = {child.tag: child for child in node}
    if len(result) != len(node):
        raise ValueError(f"Repeated child tag in {node.tag}")
    return result


def compare(original, generated, replay_tolerance=0.0):
    if not math.isfinite(replay_tolerance) or replay_tolerance < 0:
        raise ValueError("Replay tolerance must be finite and nonnegative")
    if original.tag != generated.tag:
        raise ValueError("Different UPF elements")
    ignored = {"columns", "type"}
    if original.tag == "PP_HEADER":
        ignored.update(("author", "date", "generated", "comment", "has_wfc"))
    a = {
        key: normalized(value)
        for key, value in original.attrib.items()
        if key not in ignored
    }
    b = {
        key: normalized(value)
        for key, value in generated.attrib.items()
        if key not in ignored
    }
    if a.keys() != b.keys():
        raise ValueError(f"Different attributes in {original.tag}")
    maximum = 0.0
    numeric_attributes = {"total_psenergy", "wfc_cutoff", "rho_cutoff", "pseudo_energy"}
    for key, x in a.items():
        y = b[key]
        if x == y:
            continue
        if (
            key not in numeric_attributes
            or not isinstance(x, float)
            or not isinstance(y, float)
        ):
            raise ValueError(f"Different attribute {key} in {original.tag}")
        scale = max(abs(x), abs(y))
        difference = abs(x - y) / scale
        if difference > replay_tolerance:
            raise ValueError(f"Different attribute {key} in {original.tag}")
        maximum = max(maximum, difference)
    ca, cb = children(original), children(generated)
    if original.tag == "UPF":
        for name in ("PP_INFO", "PP_FULL_WFC"):
            ca.pop(name, None)
            cb.pop(name, None)
    if ca.keys() != cb.keys():
        raise ValueError(f"Different data sections in {original.tag}")
    if ca:
        for name in ca:
            maximum = max(maximum, compare(ca[name], cb[name], replay_tolerance))
    elif original.tag != "PP_HEADER":
        va, vb = values(original), values(generated)
        if len(va) != len(vb):
            raise ValueError(f"Different array lengths in {original.tag}")
        tolerance = 0.0
        if original.tag.startswith("PP_CHI.") or original.tag == "PP_RHOATOM":
            tolerance = 64 * sys.float_info.epsilon * max(map(abs, va), default=0.0)
        scale = max(max(map(abs, va), default=0.0), max(map(abs, vb), default=0.0))
        if original.tag not in ("PP_R", "PP_RAB"):
            tolerance = max(tolerance, replay_tolerance * scale)
        difference = max((abs(x - y) for x, y in zip(va, vb)), default=0.0)
        if original.tag.startswith("PP_CHI."):
            # Atomic guess phases do not change the original USPP operators.
            difference = min(
                difference, max((abs(x + y) for x, y in zip(va, vb)), default=0.0)
            )
        if difference > tolerance:
            raise ValueError(f"Different numeric data in {original.tag}")
        if scale:
            maximum = max(maximum, difference / scale)
    return maximum


def top_level_spans(data):
    """Locate actual XML elements, preserving comments and CDATA in the output."""
    parser = expat.ParserCreate()
    depth = 0
    start = 0
    spans = {}
    root_end = None

    def opening(name, attributes):
        nonlocal depth, start
        if depth == 1:
            start = parser.CurrentByteIndex
        depth += 1

    def closing(name):
        nonlocal depth, root_end
        depth -= 1
        end = parser.CurrentByteIndex
        if depth == 1:
            if data[end : end + 2] == b"</":
                end = data.index(b">", end) + 1
            if name in spans:
                raise ValueError(f"Repeated top-level tag {name}")
            spans[name] = (start, end)
        elif depth == 0:
            root_end = end

    parser.StartElementHandler = opening
    parser.EndElementHandler = closing
    parser.Parse(data, True)
    if root_end is None:
        raise ValueError("Missing UPF closing tag")
    return spans, root_end
