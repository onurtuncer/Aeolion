// Mesh.h
//
// Minimal triangle-mesh I/O: reads/writes STL (ASCII + binary), with no
// external dependencies. This is the common substrate that both the
// direct-STL workflow and the STEP-converted workflow feed into.
//
// A "PartMesh" is meant to represent ONE aircraft component (e.g. just the
// wing, just the horizontal tail). The recommended workflow is to export
// each lifting surface (and the fuselage, if you want it along for
// bookkeeping even though it isn't lattice-meshed) as its own STL file from
// your CAD tool -- that sidesteps automatic part segmentation, which is a
// much harder and less reliable problem than slicing a single known part.

#pragma once
#include <string>
#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <limits>
#include <cstring>
#include <cstdint>
#include "Aeolion/Vlm.h"

namespace Aeolion { namespace MeshIO {

using Vlm::Vec3;

struct Triangle {
    Vec3 v0, v1, v2;
    Vec3 Normal; // as-read face normal (not guaranteed consistent/outward)
};

struct PartMesh {
    std::string Name;
    std::vector<Triangle> Tris;

    void Bbox(Vec3& lo, Vec3& hi) const {
        lo = Vec3(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
        hi = Vec3(-std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(), -std::numeric_limits<double>::max());
        auto grow = [&](const Vec3& p) {
            lo.x = std::min(lo.x, p.x); lo.y = std::min(lo.y, p.y); lo.z = std::min(lo.z, p.z);
            hi.x = std::max(hi.x, p.x); hi.y = std::max(hi.y, p.y); hi.z = std::max(hi.z, p.z);
        };
        for (const auto& t : Tris) { grow(t.v0); grow(t.v1); grow(t.v2); }
    }
};

// ------------------------------------------------------------ STL read ---

inline bool LooksLikeAsciiStl(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open " + path);
    char header[6] = {0};
    f.read(header, 5);
    std::string h(header);
    // ASCII STL files start with "solid"; binary ones usually don't, but
    // some binary files begin with "solid" too (bad exporters), so also
    // sanity-check the declared binary triangle count against file size.
    if (h.substr(0, 5) != "solid") return false;

    f.seekg(0, std::ios::end);
    std::streamsize size = f.tellg();
    f.seekg(80, std::ios::beg);
    uint32_t triCount = 0;
    f.read(reinterpret_cast<char*>(&triCount), 4);
    std::streamsize expectedBinSize = 84 + static_cast<std::streamsize>(triCount) * 50;
    if (size == expectedBinSize) return false; // actually binary despite "solid" header
    return true;
}

// Reads one or more named solids from an ASCII STL (handles multi-solid
// files: "solid A ... endsolid A  solid B ... endsolid B").
inline std::vector<PartMesh> ReadStlAscii(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open " + path);

    std::vector<PartMesh> parts;
    std::string line, tok;
    PartMesh* cur = nullptr;
    Vec3 normal;
    std::vector<Vec3> verts;

    while (std::getline(f, line)) {
        std::istringstream iss(line);
        iss >> tok;
        if (tok == "solid") {
            parts.emplace_back();
            cur = &parts.back();
            std::string nm; std::getline(iss, nm);
            // trim
            size_t a = nm.find_first_not_of(" \t\r");
            cur->Name = (a == std::string::npos) ? ("part_" + std::to_string(parts.size() - 1)) : nm.substr(a);
        } else if (tok == "facet") {
            std::string norm_kw; iss >> norm_kw >> normal.x >> normal.y >> normal.z;
            verts.clear();
        } else if (tok == "vertex") {
            Vec3 v; iss >> v.x >> v.y >> v.z;
            verts.push_back(v);
        } else if (tok == "endfacet") {
            if (verts.size() == 3 && cur) {
                cur->Tris.push_back({verts[0], verts[1], verts[2], normal});
            }
        }
    }
    if (parts.empty()) throw std::runtime_error("no 'solid' blocks found in " + path);
    return parts;
}

inline std::vector<PartMesh> ReadStlBinary(const std::string& path, const std::string& fallbackName) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open " + path);
    char header[80];
    f.read(header, 80);
    uint32_t triCount = 0;
    f.read(reinterpret_cast<char*>(&triCount), 4);

    PartMesh part;
    part.Name = fallbackName;
    part.Tris.reserve(triCount);

    for (uint32_t i = 0; i < triCount; ++i) {
        float nx, ny, nz;
        float v[9];
        f.read(reinterpret_cast<char*>(&nx), 4);
        f.read(reinterpret_cast<char*>(&ny), 4);
        f.read(reinterpret_cast<char*>(&nz), 4);
        f.read(reinterpret_cast<char*>(v), 9 * 4);
        uint16_t attrByteCount = 0;
        f.read(reinterpret_cast<char*>(&attrByteCount), 2);
        if (!f) throw std::runtime_error("unexpected EOF reading binary STL: " + path);

        Triangle t;
        t.Normal = Vec3(nx, ny, nz);
        t.v0 = Vec3(v[0], v[1], v[2]);
        t.v1 = Vec3(v[3], v[4], v[5]);
        t.v2 = Vec3(v[6], v[7], v[8]);
        part.Tris.push_back(t);
    }
    std::vector<PartMesh> parts;
    parts.push_back(std::move(part));
    return parts;
}

// name hint is used as the part name for binary STL (which has no name
// field) -- pass e.g. the filename stem.
inline std::vector<PartMesh> ReadStl(const std::string& path, const std::string& nameHint) {
    if (LooksLikeAsciiStl(path)) return ReadStlAscii(path);
    return ReadStlBinary(path, nameHint);
}

// ----------------------------------------------------------- STL write ---

inline void WriteStlBinary(const std::string& path, const PartMesh& part) {
    std::ofstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open for write: " + path);
    char header[80] = {0};
    std::string h = "generated by Mesh.h: " + part.Name;
    std::memcpy(header, h.c_str(), std::min(h.size(), size_t(79)));
    f.write(header, 80);
    uint32_t n = static_cast<uint32_t>(part.Tris.size());
    f.write(reinterpret_cast<const char*>(&n), 4);
    for (const auto& t : part.Tris) {
        float buf[12] = {
            (float)t.Normal.x, (float)t.Normal.y, (float)t.Normal.z,
            (float)t.v0.x, (float)t.v0.y, (float)t.v0.z,
            (float)t.v1.x, (float)t.v1.y, (float)t.v1.z,
            (float)t.v2.x, (float)t.v2.y, (float)t.v2.z
        };
        f.write(reinterpret_cast<const char*>(buf), sizeof(buf));
        uint16_t attr = 0;
        f.write(reinterpret_cast<const char*>(&attr), 2);
    }
}

}} // namespace Aeolion::MeshIO
