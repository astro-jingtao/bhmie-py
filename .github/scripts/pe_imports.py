"""Print the DLL import table of a Windows PE file.

Dependency-free CI diagnostic for the Windows mingw-gfortran corner job:
when an extension module fails to import with Windows' generic
"The specified module could not be found" (which names no DLL), this
script lists exactly which DLLs the module's PE import table references.
Git Bash ships no ldd and no objdump, hence this pure-stdlib parser.

Usage: python pe_imports.py FILE [FILE ...]
"""

import struct
import sys


def _rva_to_offset(rva, sections):
    for virt_addr, virt_size, raw_ptr in sections:
        if virt_addr <= rva < virt_addr + max(virt_size, 1):
            return raw_ptr + rva - virt_addr
    return None


def pe_imports(path):
    """Return the list of DLL names in the PE import table of ``path``."""
    with open(path, "rb") as fh:
        data = fh.read()

    pe_offset = struct.unpack_from("<I", data, 0x3C)[0]
    if data[pe_offset:pe_offset + 4] != b"PE\0\0":
        raise ValueError(f"{path}: not a PE file")

    num_sections = struct.unpack_from("<H", data, pe_offset + 6)[0]
    opt_header_size = struct.unpack_from("<H", data, pe_offset + 20)[0]
    opt = pe_offset + 24
    magic = struct.unpack_from("<H", data, opt)[0]

    # DataDirectory[1] (the import table) sits at opt+120 for PE32+ and
    # opt+104 for PE32 (each data-directory entry is 8 bytes).
    if magic == 0x20B:  # PE32+
        import_rva = struct.unpack_from("<I", data, opt + 120)[0]
    elif magic == 0x10B:  # PE32
        import_rva = struct.unpack_from("<I", data, opt + 104)[0]
    else:
        raise ValueError(f"{path}: unknown PE magic {magic:#x}")

    sections_start = opt + opt_header_size
    sections = []
    for i in range(num_sections):
        base = sections_start + 40 * i
        virt_size, virt_addr, _, raw_ptr = struct.unpack_from(
            "<IIII", data, base + 8
        )
        sections.append((virt_addr, virt_size, raw_ptr))

    imports = []
    desc = _rva_to_offset(import_rva, sections)
    if desc is None:
        return imports
    for i in range(1000):  # terminated by an all-zero descriptor
        entry = desc + 20 * i
        name_rva = struct.unpack_from("<I", data, entry + 12)[0]
        if name_rva == 0:
            break
        name_off = _rva_to_offset(name_rva, sections)
        end = data.index(b"\0", name_off)
        imports.append(data[name_off:end].decode("ascii"))
    return imports


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    for path in sys.argv[1:]:
        print(f"{path}:")
        for dll in pe_imports(path):
            print(f"  {dll}")
