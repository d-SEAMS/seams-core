//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#include <pore.hpp>

#include <generic.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <utility>

namespace pore {
namespace {

constexpr double kPi = std::numbers::pi_v<double>;

struct Cell {
  double lx = 0, ly = 0, lz = 0;
  double xy = 0, xz = 0, yz = 0;
  double ox = 0, oy = 0, oz = 0;
  double volume = 0;
  double length[3] = {0, 0, 0};
  double lattice[3][3] = {};
  bool ok = false;
};

struct P2 {
  double x = 0;
  double y = 0;
};

double wrap01(double s) {
  double w = s - std::floor(s);
  if (w >= 1.0) {
    w = 0.0;
  }
  if (w < 0.0) {
    w = 0.0;
  }
  return w;
}

double wrapHalf(double s) {
  double w = wrap01(s);
  if (w >= 0.5) {
    w -= 1.0;
  }
  return w;
}

double dot3(const double a[3], const double b[3]) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

void cross3(const double a[3], const double b[3], double out[3]) {
  out[0] = a[1] * b[2] - a[2] * b[1];
  out[1] = a[2] * b[0] - a[0] * b[2];
  out[2] = a[0] * b[1] - a[1] * b[0];
}

double norm3(const double a[3]) { return std::sqrt(dot3(a, a)); }

Cell cellOf(const Cloud &cloud) {
  const gen::FracBox box = gen::makeFracBox(cloud);
  Cell cell;
  if (!box.ok) {
    return cell;
  }
  cell.lx = box.lx;
  cell.ly = box.ly;
  cell.lz = box.lz;
  cell.xy = box.xy;
  cell.xz = box.xz;
  cell.yz = box.yz;
  cell.ox = box.ox;
  cell.oy = box.oy;
  cell.oz = box.oz;
  cell.volume = box.lx * box.ly * box.lz;
  cell.lattice[0][0] = box.lx;
  cell.lattice[1][0] = box.xy;
  cell.lattice[1][1] = box.ly;
  cell.lattice[2][0] = box.xz;
  cell.lattice[2][1] = box.yz;
  cell.lattice[2][2] = box.lz;
  cell.length[0] = box.lx;
  cell.length[1] = std::hypot(box.xy, box.ly);
  cell.length[2] = std::hypot(std::hypot(box.xz, box.yz), box.lz);
  cell.ok = cell.volume > 0.0 && cell.length[0] > 0.0 && cell.length[1] > 0.0 &&
            cell.length[2] > 0.0;
  return cell;
}

void requireCell(const Cell &cell) {
  if (!cell.ok) {
    throw BoxError("box is not a positive cell");
  }
}

void toFrac(const Cell &cell, double x, double y, double z, double s[3]) {
  s[2] = (z - cell.oz) / cell.lz;
  s[1] = (y - cell.oy - cell.yz * s[2]) / cell.ly;
  s[0] = (x - cell.ox - cell.xy * s[1] - cell.xz * s[2]) / cell.lx;
}

void fromFrac(const Cell &cell, const double s[3], double r[3]) {
  r[0] = cell.ox + cell.lx * s[0] + cell.xy * s[1] + cell.xz * s[2];
  r[1] = cell.oy + cell.ly * s[1] + cell.yz * s[2];
  r[2] = cell.oz + cell.lz * s[2];
}

void requireBins(int bins) {
  if (bins < 1) {
    throw std::out_of_range("bins must be positive");
  }
}

void requireOrigin(double origin) {
  if (!std::isfinite(origin)) {
    throw std::invalid_argument("origin is not finite");
  }
}

double massAt(const double *masses, int index) {
  if (masses == nullptr) {
    return 1.0;
  }
  const double mass = masses[index];
  if (!std::isfinite(mass) || mass < 0.0) {
    throw std::invalid_argument("mass is not a finite non-negative number");
  }
  return mass;
}

struct Sample {
  double s[3];
  double mass;
};

std::vector<Sample> select(const Cloud &cloud, const Cell &cell, int type,
                           const double *masses) {
  std::vector<Sample> out;
  out.reserve(cloud.pts.size());
  for (int i = 0; i < static_cast<int>(cloud.pts.size()); ++i) {
    const auto &pt = cloud.pts[static_cast<std::size_t>(i)];
    if (type != 0 && pt.type != type) {
      continue;
    }
    Sample sample;
    toFrac(cell, pt.x, pt.y, pt.z, sample.s);
    sample.s[0] = wrap01(sample.s[0]);
    sample.s[1] = wrap01(sample.s[1]);
    sample.s[2] = wrap01(sample.s[2]);
    sample.mass = massAt(masses, i);
    out.push_back(sample);
  }
  return out;
}

Profile blank(int bins) {
  Profile out;
  out.coordinate.assign(static_cast<std::size_t>(bins), 0.0);
  out.number_density.assign(static_cast<std::size_t>(bins), 0.0);
  out.mass_density.assign(static_cast<std::size_t>(bins), 0.0);
  out.volume.assign(static_cast<std::size_t>(bins), 0.0);
  return out;
}

void finish(Profile &out, const std::vector<double> &count, const std::vector<double> &mass) {
  for (std::size_t i = 0; i < out.volume.size(); ++i) {
    const double volume = out.volume[i];
    if (volume > 0.0) {
      out.number_density[i] = count[i] / volume;
      out.mass_density[i] = mass[i] / volume;
    }
  }
}

int binOf(double u, int bins) {
  int bin = static_cast<int>(std::floor(u * static_cast<double>(bins)));
  if (bin >= bins) {
    bin = bins - 1;
  }
  if (bin < 0) {
    bin = 0;
  }
  return bin;
}

// Volume of {s in [0,1]^m : g·s <= t} for g_i > 0, embedded in the unit cube.
double positiveCdf(const double *g, int m, double t) {
  if (m <= 0) {
    return t >= 0.0 ? 1.0 : 0.0;
  }
  double prod = 1.0;
  for (int i = 0; i < m; ++i) {
    prod *= g[i];
  }
  double fact = 1.0;
  for (int i = 1; i <= m; ++i) {
    fact *= static_cast<double>(i);
  }
  const int nvert = 1 << m;
  double acc = 0.0;
  for (int mask = 0; mask < nvert; ++mask) {
    int bits = 0;
    double shift = 0.0;
    for (int i = 0; i < m; ++i) {
      if ((mask & (1 << i)) != 0) {
        bits += 1;
        shift += g[i];
      }
    }
    const double arg = t - shift;
    if (arg > 0.0) {
      double term = 1.0;
      for (int p = 0; p < m; ++p) {
        term *= arg;
      }
      acc += (bits % 2 == 0) ? term : -term;
    }
  }
  double vol = acc / (fact * prod);
  if (vol < 0.0) {
    vol = 0.0;
  }
  if (vol > 1.0) {
    vol = 1.0;
  }
  return vol;
}

double cubeCdf(const double g[3], double x) {
  double shift = 0.0;
  double gp[3];
  int m = 0;
  for (int i = 0; i < 3; ++i) {
    if (g[i] < 0.0) {
      shift += g[i];
      gp[m++] = -g[i];
    } else if (g[i] > 0.0) {
      gp[m++] = g[i];
    }
  }
  return positiveCdf(gp, m, x - shift);
}

double cubeMeasure(const double g[3], double lo, double hi) {
  if (!(hi > lo)) {
    return 0.0;
  }
  return cubeCdf(g, hi) - cubeCdf(g, lo);
}

struct DirectionGeom {
  double g[3];
  double shift = 0;
  double span = 0;
};

DirectionGeom directionGeom(const Cell &cell, const double direction[3]) {
  const double scale = norm3(direction);
  if (!(scale > 0.0) || !std::isfinite(scale)) {
    throw std::invalid_argument("direction is zero");
  }
  const double u[3] = {direction[0] / scale, direction[1] / scale, direction[2] / scale};
  DirectionGeom geom;
  for (int i = 0; i < 3; ++i) {
    geom.g[i] = dot3(cell.lattice[i], u);
  }
  double positive = 0.0;
  for (int i = 0; i < 3; ++i) {
    if (geom.g[i] < 0.0) {
      geom.shift += geom.g[i];
    } else {
      positive += geom.g[i];
    }
  }
  geom.span = positive - geom.shift;
  if (!(geom.span > 0.0)) {
    throw std::invalid_argument("direction is zero");
  }
  return geom;
}

Profile project(const Cloud &cloud, int type, int bins, double origin, const double *masses,
                const DirectionGeom &geom) {
  requireBins(bins);
  requireOrigin(origin);
  const Cell cell = cellOf(cloud);
  requireCell(cell);
  const double org = wrap01(origin);
  Profile out = blank(bins);
  std::vector<double> count(static_cast<std::size_t>(bins), 0.0);
  std::vector<double> mass(static_cast<std::size_t>(bins), 0.0);
  for (int i = 0; i < bins; ++i) {
    const double a = org + static_cast<double>(i) / static_cast<double>(bins);
    const double b = org + static_cast<double>(i + 1) / static_cast<double>(bins);
    double measure = 0.0;
    // A bin that starts at or past the wrap lives entirely in the next
    // cell image. Measuring from 0 up to b-1 counts every earlier slab
    // again, so the volumes no longer sum to the cell.
    if (a >= 1.0) {
      measure = cubeMeasure(geom.g, geom.shift + (a - 1.0) * geom.span,
                             geom.shift + (b - 1.0) * geom.span);
    } else if (b <= 1.0) {
      measure = cubeMeasure(geom.g, geom.shift + a * geom.span, geom.shift + b * geom.span);
    } else {
      measure = cubeMeasure(geom.g, geom.shift + a * geom.span, geom.shift + geom.span);
      measure += cubeMeasure(geom.g, geom.shift, geom.shift + (b - 1.0) * geom.span);
    }
    out.volume[static_cast<std::size_t>(i)] = measure * cell.volume;
    out.coordinate[static_cast<std::size_t>(i)] =
        (static_cast<double>(i) + 0.5) / static_cast<double>(bins);
  }
  for (const Sample &sample : select(cloud, cell, type, masses)) {
    const double xi = geom.g[0] * sample.s[0] + geom.g[1] * sample.s[1] + geom.g[2] * sample.s[2];
    const double t = (xi - geom.shift) / geom.span;
    const int bin = binOf(wrap01(t - org), bins);
    count[static_cast<std::size_t>(bin)] += 1.0;
    mass[static_cast<std::size_t>(bin)] += sample.mass;
  }
  finish(out, count, mass);
  return out;
}

struct Plane {
  int axis = 2;
  int iu = 0;
  int iv = 1;
  double e1[3] = {};
  double e2[3] = {};
  double eaxis[3] = {};
  P2 eu;
  P2 ev;
  double axial = 0;
  double face = 0;
  double far = 0;
  double cover = 0;
  std::vector<P2> polygon;
  std::vector<P2> voronoi;
};

std::vector<P2> voronoiCell(const Plane &plane);

int parallelAxis(const Cell &cell, const double direction[3]) {
  const double dw = direction[2] / cell.lz;
  const double dv = (direction[1] - cell.yz * dw) / cell.ly;
  const double du = (direction[0] - cell.xy * dv - cell.xz * dw) / cell.lx;
  const double frac[3] = {du, dv, dw};
  const double scale =
      std::max(std::abs(frac[0]), std::max(std::abs(frac[1]), std::abs(frac[2])));
  if (!(scale > 0.0) || !std::isfinite(scale)) {
    return -1;
  }
  int nonzero = -1;
  int count = 0;
  for (int i = 0; i < 3; ++i) {
    if (std::abs(frac[i]) > 1e-12 * scale) {
      nonzero = i;
      count += 1;
    }
  }
  return count == 1 ? nonzero : -1;
}

Plane planeOf(const Cell &cell, int axis) {
  if (axis < 0 || axis > 2) {
    throw std::out_of_range("axis is 0, 1, or 2");
  }
  Plane plane;
  plane.axis = axis;
  plane.iu = (axis + 1) % 3;
  plane.iv = (axis + 2) % 3;
  const double *along = cell.lattice[axis];
  const double axial2 = dot3(along, along);
  plane.axial = std::sqrt(axial2);
  for (int k = 0; k < 3; ++k) {
    plane.eaxis[k] = along[k] / plane.axial;
  }
  double raw_u[3];
  double raw_v[3];
  const double pu = dot3(cell.lattice[plane.iu], along) / axial2;
  const double pv = dot3(cell.lattice[plane.iv], along) / axial2;
  for (int k = 0; k < 3; ++k) {
    raw_u[k] = cell.lattice[plane.iu][k] - pu * along[k];
    raw_v[k] = cell.lattice[plane.iv][k] - pv * along[k];
  }
  const double nu = norm3(raw_u);
  const double nv = norm3(raw_v);
  if (!(nu > 0.0) || !(nv > 0.0)) {
    throw BoxError("box is not a positive cell");
  }
  for (int k = 0; k < 3; ++k) {
    plane.e1[k] = raw_u[k] / nu;
  }
  cross3(plane.eaxis, plane.e1, plane.e2);
  const double n2 = norm3(plane.e2);
  if (!(n2 > 0.0)) {
    throw BoxError("box is not a positive cell");
  }
  for (int k = 0; k < 3; ++k) {
    plane.e2[k] /= n2;
  }
  plane.eu = {dot3(raw_u, plane.e1), dot3(raw_u, plane.e2)};
  plane.ev = {dot3(raw_v, plane.e1), dot3(raw_v, plane.e2)};
  const double signed_area = plane.eu.x * plane.ev.y - plane.eu.y * plane.ev.x;
  plane.face = std::abs(signed_area);
  if (!(plane.face > 0.0)) {
    throw BoxError("box is not a positive cell");
  }
  plane.polygon = {
      {-0.5 * plane.eu.x - 0.5 * plane.ev.x, -0.5 * plane.eu.y - 0.5 * plane.ev.y},
      {0.5 * plane.eu.x - 0.5 * plane.ev.x, 0.5 * plane.eu.y - 0.5 * plane.ev.y},
      {0.5 * plane.eu.x + 0.5 * plane.ev.x, 0.5 * plane.eu.y + 0.5 * plane.ev.y},
      {-0.5 * plane.eu.x + 0.5 * plane.ev.x, -0.5 * plane.eu.y + 0.5 * plane.ev.y},
  };
  if (signed_area < 0.0) {
    std::swap(plane.polygon[1], plane.polygon[3]);
  }
  plane.far = 0.0;
  for (const P2 &p : plane.polygon) {
    plane.far = std::max(plane.far, std::hypot(p.x, p.y));
  }
  plane.voronoi = voronoiCell(plane);
  plane.cover = 0.0;
  for (const P2 &p : plane.voronoi) {
    plane.cover = std::max(plane.cover, std::hypot(p.x, p.y));
  }
  if (!(plane.cover > 0.0)) {
    throw BoxError("box is not a positive cell");
  }
  return plane;
}

double cross2(P2 a, P2 b) { return a.x * b.y - a.y * b.x; }

double dot2(P2 a, P2 b) { return a.x * b.x + a.y * b.y; }

double shoelace(const std::vector<P2> &poly) {
  double acc = 0.0;
  const int n = static_cast<int>(poly.size());
  for (int i = 0; i < n; ++i) {
    const P2 &a = poly[static_cast<std::size_t>(i)];
    const P2 &b = poly[static_cast<std::size_t>((i + 1) % n)];
    acc += cross2(a, b);
  }
  return 0.5 * acc;
}

// Signed area of a convex polygon intersected with the disk at the origin.
// The disk is convex, so an edge with both ends inside lies inside the disk.
double diskArea(const std::vector<P2> &poly, double radius, double far, double face) {
  if (!(radius > 0.0)) {
    return 0.0;
  }
  if (radius >= far) {
    return face;
  }
  const double r2 = radius * radius;
  const int n = static_cast<int>(poly.size());
  double area = 0.0;
  for (int i = 0; i < n; ++i) {
    const P2 a = poly[static_cast<std::size_t>(i)];
    const P2 b = poly[static_cast<std::size_t>((i + 1) % n)];
    const bool ain = dot2(a, a) <= r2;
    const bool bin = dot2(b, b) <= r2;
    P2 delta{b.x - a.x, b.y - a.y};
    const double aa = dot2(delta, delta);
    const double bb = 2.0 * dot2(a, delta);
    const double cc = dot2(a, a) - r2;
    std::vector<P2> hits;
    const double disc = bb * bb - 4.0 * aa * cc;
    if (aa > 0.0 && disc >= 0.0) {
      const double root = std::sqrt(disc);
      const double q = -0.5 * (bb + std::copysign(root, bb));
      double roots[2];
      int nroot = 0;
      if (q != 0.0) {
        roots[nroot++] = q / aa;
        roots[nroot++] = cc / q;
      } else {
        roots[nroot++] = root == 0.0 ? 0.0 : -bb / (2.0 * aa);
      }
      std::sort(roots, roots + nroot);
      // A root within 1e-12 of an endpoint is the crossing at that vertex.
      // Dropping it makes the edge a sector, and a Voronoi cell whose
      // vertices share one radius then sums to the full disk. Clamp the
      // root onto the segment and keep the point.
      bool have = false;
      double accepted = 0.0;
      for (int k = 0; k < nroot; ++k) {
        const double along = roots[k];
        if (along < -1e-12 || along > 1.0 + 1e-12) {
          continue;
        }
        const double t = std::clamp(along, 0.0, 1.0);
        if (have && t - accepted < 1e-15) {
          continue;
        }
        have = true;
        accepted = t;
        hits.push_back(P2{a.x + t * delta.x, a.y + t * delta.y});
      }
    }
    auto tri = [](P2 u, P2 v) { return 0.5 * cross2(u, v); };
    auto sector = [&](P2 u, P2 v) {
      return 0.5 * r2 * std::atan2(cross2(u, v), dot2(u, v));
    };
    if (ain && bin) {
      area += tri(a, b);
    } else if (ain) {
      // a is inside and b is outside. hits.back() is the exit. hits.front()
      // can be a itself once an endpoint root is kept, and using it would
      // throw the exit away.
      const P2 h = hits.empty() ? b : hits.back();
      area += tri(a, h) + sector(h, b);
    } else if (bin) {
      const P2 h = hits.empty() ? a : hits.front();
      area += sector(a, h) + tri(h, b);
    } else if (hits.size() >= 2) {
      area += sector(a, hits.front()) + tri(hits.front(), hits.back()) + sector(hits.back(), b);
    } else {
      // The edge misses the disk, or it touches at one interior point.
      // Either way the chord has no area.
      area += sector(a, b);
    }
  }
  return area;
}

bool leftOf(P2 p, P2 a, P2 b) {
  return cross2(P2{b.x - a.x, b.y - a.y}, P2{p.x - a.x, p.y - a.y}) >= -1e-14;
}

P2 lineHit(P2 p, P2 q, P2 a, P2 b) {
  const double rx = q.x - p.x;
  const double ry = q.y - p.y;
  const double sx = b.x - a.x;
  const double sy = b.y - a.y;
  const double den = rx * sy - ry * sx;
  const double t = ((a.x - p.x) * sy - (a.y - p.y) * sx) / den;
  return {p.x + t * rx, p.y + t * ry};
}

std::vector<P2> clipConvex(std::vector<P2> subject, const std::vector<P2> &clip) {
  const int nc = static_cast<int>(clip.size());
  for (int e = 0; e < nc; ++e) {
    if (subject.empty()) {
      break;
    }
    const P2 a = clip[static_cast<std::size_t>(e)];
    const P2 b = clip[static_cast<std::size_t>((e + 1) % nc)];
    std::vector<P2> next;
    const int ns = static_cast<int>(subject.size());
    for (int i = 0; i < ns; ++i) {
      const P2 cur = subject[static_cast<std::size_t>(i)];
      const P2 prev = subject[static_cast<std::size_t>((i + ns - 1) % ns)];
      const bool cin = leftOf(cur, a, b);
      const bool pin = leftOf(prev, a, b);
      const double rx = cur.x - prev.x;
      const double ry = cur.y - prev.y;
      const double sx = b.x - a.x;
      const double sy = b.y - a.y;
      const double den = rx * sy - ry * sx;
      const double scale = std::hypot(rx, ry) * std::hypot(sx, sy);
      const bool crosses = std::abs(den) > 1e-14 * std::max(scale, 1.0);
      if (cin && pin) {
        next.push_back(cur);
      } else if (cin && !pin) {
        if (crosses) {
          next.push_back(lineHit(prev, cur, a, b));
        }
        next.push_back(cur);
      } else if (!cin && pin && crosses) {
        next.push_back(lineHit(prev, cur, a, b));
      }
    }
    subject.swap(next);
  }
  return subject;
}

// Voronoi cell of the perpendicular lattice. Its points are closer to
// the pore axis than to any periodic image of that axis, and its area
// is the face area. Radial shells are disks clipped to this cell, so
// the shell volume matches the nearest-image radius.
std::vector<P2> voronoiCell(const Plane &plane) {
  const double box = plane.far * 4.0 + 1.0;
  std::vector<P2> poly = {{-box, -box}, {box, -box}, {box, box}, {-box, box}};
  const double limit = 2.0 * plane.far + 1.0;
  const double lu = std::hypot(plane.eu.x, plane.eu.y);
  const double lv = std::hypot(plane.ev.x, plane.ev.y);
  const double span = std::max(plane.face, 1e-30);
  const int ni = static_cast<int>(std::ceil(limit * lv / span)) + 2;
  const int nj = static_cast<int>(std::ceil(limit * lu / span)) + 2;
  struct Site {
    double r2;
    int i;
    int j;
  };
  std::vector<Site> sites;
  for (int i = -ni; i <= ni; ++i) {
    for (int j = -nj; j <= nj; ++j) {
      if (i == 0 && j == 0) {
        continue;
      }
      const P2 w{static_cast<double>(i) * plane.eu.x + static_cast<double>(j) * plane.ev.x,
                 static_cast<double>(i) * plane.eu.y + static_cast<double>(j) * plane.ev.y};
      const double r2 = dot2(w, w);
      if (r2 > 0.0 && r2 <= limit * limit) {
        sites.push_back(Site{r2, i, j});
      }
    }
  }
  std::sort(sites.begin(), sites.end(), [](const Site &a, const Site &b) {
    if (a.r2 != b.r2) {
      return a.r2 < b.r2;
    }
    if (a.i != b.i) {
      return a.i < b.i;
    }
    return a.j < b.j;
  });
  for (const Site &site : sites) {
    const P2 w{static_cast<double>(site.i) * plane.eu.x + static_cast<double>(site.j) * plane.ev.x,
               static_cast<double>(site.i) * plane.eu.y + static_cast<double>(site.j) * plane.ev.y};
    const P2 mid{0.5 * w.x, 0.5 * w.y};
    // rot90ccw(w) = (-wy, wx). The interior (the origin side) is to the left.
    P2 dir{-w.y, w.x};
    const double dir_n = std::hypot(dir.x, dir.y);
    if (!(dir_n > 0.0)) {
      continue;
    }
    dir.x /= dir_n;
    dir.y /= dir_n;
    const P2 inward{-dir.y, dir.x};
    const double reach = box * 8.0 + 1.0;
    const P2 a{mid.x - reach * dir.x, mid.y - reach * dir.y};
    const P2 b{mid.x + reach * dir.x, mid.y + reach * dir.y};
    const P2 c{b.x + reach * inward.x, b.y + reach * inward.y};
    const P2 d{a.x + reach * inward.x, a.y + reach * inward.y};
    poly = clipConvex(poly, std::vector<P2>{a, b, c, d});
    if (poly.empty()) {
      break;
    }
  }
  return poly;
}

double wedgeArea(const Plane &plane, double theta0, double theta1) {
  const double span = theta1 - theta0;
  if (!(span > 0.0)) {
    return 0.0;
  }
  if (span >= 2.0 * kPi - 1e-12) {
    return plane.face;
  }
  auto piece = [&](double a0, double a1) {
    const double reach = plane.far * 4.0 + 1.0;
    std::vector<P2> wedge = {
        {0.0, 0.0},
        {reach * std::cos(a0), reach * std::sin(a0)},
        {reach * std::cos(a1), reach * std::sin(a1)},
    };
    return std::abs(shoelace(clipConvex(plane.polygon, wedge)));
  };
  // The clipping triangle covers the sector only while its chord stays
  // outside the face. reach is about 4*far, so that fails above about
  // 151 degrees, and at span = pi the three points are collinear and
  // the triangle has no area. Pieces of at most pi/2 keep the chord
  // outside the face.
  constexpr double kPiece = 0.5 * kPi;
  if (span <= kPiece) {
    return piece(theta0, theta1);
  }
  const int pieces = static_cast<int>(std::ceil(span / kPiece - 1e-15));
  double area = 0.0;
  for (int i = 0; i < pieces; ++i) {
    const double a0 = theta0 + span * (static_cast<double>(i) / static_cast<double>(pieces));
    const double a1 = theta0 + span * (static_cast<double>(i + 1) / static_cast<double>(pieces));
    area += piece(a0, a1);
  }
  return area;
}

// Position in the centred perpendicular face. Each face coordinate is
// wrapped into [-1/2, 1/2). The azimuth and the automatic centre use
// this map. The radial coordinate does not: past the inscribed radius
// it follows the nearest image of the pore axis.
void perpOf(const Plane &plane, const double s[3], double cu, double cv, double &x, double &y) {
  const double du = wrapHalf(s[plane.iu] - cu);
  const double dv = wrapHalf(s[plane.iv] - cv);
  x = du * plane.eu.x + dv * plane.ev.x;
  y = du * plane.eu.y + dv * plane.ev.y;
}

void imageSearch(const Plane &plane, int &ni, int &nj) {
  const double limit = 2.0 * plane.far + 1.0;
  const double lu = std::hypot(plane.eu.x, plane.eu.y);
  const double lv = std::hypot(plane.ev.x, plane.ev.y);
  const double span = std::max(plane.face, 1e-30);
  ni = static_cast<int>(std::ceil(limit * lv / span)) + 2;
  nj = static_cast<int>(std::ceil(limit * lu / span)) + 2;
}

// Displacement from the nearest periodic image of the pore axis.
// (du, dv) are the centred fractional coordinates, already in [-1/2, 1/2).
void nearestAxis(const Plane &plane, double du, double dv, double &x, double &y) {
  int ni = 0;
  int nj = 0;
  imageSearch(plane, ni, nj);
  double best = std::numeric_limits<double>::infinity();
  double bx = 0.0;
  double by = 0.0;
  for (int i = -ni; i <= ni; ++i) {
    for (int j = -nj; j <= nj; ++j) {
      const double px = (du + static_cast<double>(i)) * plane.eu.x +
                        (dv + static_cast<double>(j)) * plane.ev.x;
      const double py = (du + static_cast<double>(i)) * plane.eu.y +
                        (dv + static_cast<double>(j)) * plane.ev.y;
      const double r2 = px * px + py * py;
      if (r2 < best) {
        best = r2;
        bx = px;
        by = py;
      }
    }
  }
  x = bx;
  y = by;
}

struct Centre2 {
  double u;
  double v;
};

double meanRadius2(const std::vector<Sample> &walls, const Plane &plane, double u, double v) {
  double acc = 0.0;
  for (const Sample &wall : walls) {
    double x = 0.0;
    double y = 0.0;
    perpOf(plane, wall.s, u, v, x, y);
    acc += x * x + y * y;
  }
  return walls.empty() ? 0.0 : acc / static_cast<double>(walls.size());
}

Centre2 iterateCentre(const std::vector<Sample> &walls, const Plane &plane, double u, double v) {
  const double det = plane.eu.x * plane.ev.y - plane.eu.y * plane.ev.x;
  for (int step = 0; step < 40; ++step) {
    double mx = 0.0;
    double my = 0.0;
    for (const Sample &wall : walls) {
      double x = 0.0;
      double y = 0.0;
      perpOf(plane, wall.s, u, v, x, y);
      mx += x;
      my += y;
    }
    const double inv = 1.0 / static_cast<double>(walls.size());
    mx *= inv;
    my *= inv;
    const double du = (mx * plane.ev.y - my * plane.ev.x) / det;
    const double dv = (plane.eu.x * my - plane.eu.y * mx) / det;
    u = wrap01(u + du);
    v = wrap01(v + dv);
    if (mx * mx + my * my <= 1e-24 * plane.far * plane.far) {
      break;
    }
  }
  return {u, v};
}

double circularMean(const std::vector<Sample> &walls, int axis) {
  double sx = 0.0;
  double sy = 0.0;
  for (const Sample &wall : walls) {
    const double ang = 2.0 * kPi * wall.s[axis];
    sx += std::cos(ang);
    sy += std::sin(ang);
  }
  const double inv = 1.0 / static_cast<double>(walls.size());
  return wrap01(std::atan2(sy * inv, sx * inv) / (2.0 * kPi));
}

Centre2 autoCentre(const std::vector<Sample> &walls, const Plane &plane) {
  if (walls.empty()) {
    throw std::invalid_argument("no framework atoms");
  }
  const Centre2 seeds[5] = {
      {0.0, 0.0},
      {0.5, 0.5},
      {0.0, 0.5},
      {0.5, 0.0},
      {circularMean(walls, plane.iu), circularMean(walls, plane.iv)},
  };
  Centre2 best = seeds[0];
  double best_score = 0.0;
  bool have = false;
  for (const Centre2 &seed : seeds) {
    const Centre2 fit = iterateCentre(walls, plane, seed.u, seed.v);
    const double score = meanRadius2(walls, plane, fit.u, fit.v);
    const bool closer = !have || score < best_score - 1e-12 * plane.far * plane.far ||
                        (std::abs(score - best_score) <= 1e-12 * plane.far * plane.far &&
                         (fit.u < best.u || (fit.u == best.u && fit.v < best.v)));
    if (closer) {
      best = fit;
      best_score = score;
      have = true;
    }
  }
  return best;
}

int autoAxis(const Cell &cell) {
  int best = 0;
  for (int axis = 1; axis < 3; ++axis) {
    if (cell.length[axis] > cell.length[best]) {
      best = axis;
    }
  }
  return best;
}

int resolveAxis(const Cell &cell, int axis, const double *direction) {
  if (direction != nullptr) {
    const int parallel = parallelAxis(cell, direction);
    if (parallel < 0) {
      throw std::invalid_argument("pore axis must be a cell edge");
    }
    return parallel;
  }
  if (axis == -1) {
    return autoAxis(cell);
  }
  if (axis < 0 || axis > 2) {
    throw std::out_of_range("axis is -1, 0, 1, or 2");
  }
  return axis;
}

void writeCentre(const Cell &cell, const Plane &plane, double u, double v, double w, Frame &frame) {
  frame.centre_frac[plane.iu] = wrap01(u);
  frame.centre_frac[plane.iv] = wrap01(v);
  frame.centre_frac[plane.axis] = wrap01(w);
  fromFrac(cell, frame.centre_frac, frame.centre);
  for (int k = 0; k < 3; ++k) {
    frame.axis_unit[k] = plane.eaxis[k];
  }
  frame.axis_index = plane.axis;
}

double angleOf(double y, double x) {
  double theta = std::atan2(y, x);
  if (theta >= kPi) {
    theta = -kPi;
  }
  return theta;
}

} // namespace

double clippedDiskArea(const std::vector<std::array<double, 2>> &polygon, double radius,
                       double cover, double full) {
  std::vector<P2> poly;
  poly.reserve(polygon.size());
  for (const auto &point : polygon) {
    poly.push_back(P2{point[0], point[1]});
  }
  return diskArea(poly, radius, cover, full);
}

Profile fractionalAxis(const Cloud &cloud, int type, int bins, int axis, double origin,
                       const double *masses) {
  requireBins(bins);
  requireOrigin(origin);
  if (axis < 0 || axis > 2) {
    throw std::out_of_range("axis is 0, 1, or 2");
  }
  const Cell cell = cellOf(cloud);
  requireCell(cell);
  const double org = wrap01(origin);
  Profile out = blank(bins);
  const double slab = cell.volume / static_cast<double>(bins);
  std::vector<double> count(static_cast<std::size_t>(bins), 0.0);
  std::vector<double> mass(static_cast<std::size_t>(bins), 0.0);
  for (int i = 0; i < bins; ++i) {
    out.volume[static_cast<std::size_t>(i)] = slab;
    out.coordinate[static_cast<std::size_t>(i)] =
        (static_cast<double>(i) + 0.5) / static_cast<double>(bins);
  }
  for (const Sample &sample : select(cloud, cell, type, masses)) {
    const int bin = binOf(wrap01(sample.s[axis] - org), bins);
    count[static_cast<std::size_t>(bin)] += 1.0;
    mass[static_cast<std::size_t>(bin)] += sample.mass;
  }
  finish(out, count, mass);
  return out;
}

Profile alongDirection(const Cloud &cloud, int type, int bins, const double direction[3],
                       double origin, const double *masses) {
  if (direction == nullptr) {
    throw std::invalid_argument("direction is zero");
  }
  const Cell cell = cellOf(cloud);
  requireCell(cell);
  return project(cloud, type, bins, origin, masses, directionGeom(cell, direction));
}

Frame poreFrame(const Cloud &cloud, int framework_type, int axis, const double *direction,
                const double *centre_frac) {
  const Cell cell = cellOf(cloud);
  requireCell(cell);
  const int resolved = resolveAxis(cell, axis, direction);
  const Plane plane = planeOf(cell, resolved);
  const auto walls = select(cloud, cell, framework_type, nullptr);
  double u = 0.0;
  double v = 0.0;
  double w = 0.0;
  if (centre_frac != nullptr) {
    u = centre_frac[plane.iu];
    v = centre_frac[plane.iv];
    w = centre_frac[plane.axis];
    if (!std::isfinite(u) || !std::isfinite(v) || !std::isfinite(w)) {
      throw std::invalid_argument("centre is not finite");
    }
  } else {
    const Centre2 fit = autoCentre(walls, plane);
    u = fit.u;
    v = fit.v;
    w = circularMean(walls, plane.axis);
  }
  Frame frame;
  writeCentre(cell, plane, u, v, w, frame);
  return frame;
}

Profiles poreProfiles(const Cloud &cloud, int type, int radial_bins, int azimuthal_bins,
                      int axial_bins, double r_max, int axis, const double *direction,
                      const double *centre_frac, int framework_type, double origin,
                      const double *masses) {
  requireBins(radial_bins);
  requireBins(azimuthal_bins);
  requireBins(axial_bins);
  requireOrigin(origin);
  if (!std::isfinite(r_max) || r_max < 0.0) {
    throw std::out_of_range("r_max must be non-negative");
  }
  const Cell cell = cellOf(cloud);
  requireCell(cell);
  const int resolved = resolveAxis(cell, axis, direction);
  const Plane plane = planeOf(cell, resolved);
  double cu = 0.0;
  double cv = 0.0;
  if (centre_frac != nullptr) {
    if (!std::isfinite(centre_frac[0]) || !std::isfinite(centre_frac[1]) ||
        !std::isfinite(centre_frac[2])) {
      throw std::invalid_argument("centre is not finite");
    }
    cu = centre_frac[plane.iu];
    cv = centre_frac[plane.iv];
  } else {
    const Frame frame = poreFrame(cloud, framework_type, resolved, nullptr, nullptr);
    cu = frame.centre_frac[plane.iu];
    cv = frame.centre_frac[plane.iv];
  }
  const double reach = r_max > 0.0 ? r_max : plane.far;
  if (!(reach > 0.0)) {
    throw BoxError("box is not a positive cell");
  }
  Profiles out;
  out.radial = blank(radial_bins);
  out.azimuthal = blank(azimuthal_bins);
  out.axial = fractionalAxis(cloud, type, axial_bins, resolved, origin, masses);
  const double dr = reach / static_cast<double>(radial_bins);
  const double dtheta = 2.0 * kPi / static_cast<double>(azimuthal_bins);
  std::vector<double> rcount(static_cast<std::size_t>(radial_bins), 0.0);
  std::vector<double> rmass(static_cast<std::size_t>(radial_bins), 0.0);
  std::vector<double> tcount(static_cast<std::size_t>(azimuthal_bins), 0.0);
  std::vector<double> tmass(static_cast<std::size_t>(azimuthal_bins), 0.0);
  for (int i = 0; i < radial_bins; ++i) {
    const double inner = dr * static_cast<double>(i);
    const double outer = dr * static_cast<double>(i + 1);
    // plane.far is the centred-face vertex radius, which is at least the
    // Voronoi covering radius. The early return then hands back the exact
    // face area once the disk covers the centred face.
    double area = diskArea(plane.voronoi, outer, plane.far, plane.face) -
                  diskArea(plane.voronoi, inner, plane.far, plane.face);
    if (area < 0.0) {
      area = 0.0;
    }
    out.radial.volume[static_cast<std::size_t>(i)] = area * plane.axial;
    out.radial.coordinate[static_cast<std::size_t>(i)] = 0.5 * (inner + outer);
  }
  for (int i = 0; i < azimuthal_bins; ++i) {
    const double a0 = -kPi + dtheta * static_cast<double>(i);
    const double a1 = -kPi + dtheta * static_cast<double>(i + 1);
    out.azimuthal.volume[static_cast<std::size_t>(i)] = wedgeArea(plane, a0, a1) * plane.axial;
    out.azimuthal.coordinate[static_cast<std::size_t>(i)] = 0.5 * (a0 + a1);
  }
  for (const Sample &sample : select(cloud, cell, type, masses)) {
    double x = 0.0;
    double y = 0.0;
    perpOf(plane, sample.s, cu, cv, x, y);
    const double du = wrapHalf(sample.s[plane.iu] - cu);
    const double dv = wrapHalf(sample.s[plane.iv] - cv);
    double xr = 0.0;
    double yr = 0.0;
    nearestAxis(plane, du, dv, xr, yr);
    const double radius = std::hypot(xr, yr);
    if (radius <= reach) {
      int bin = static_cast<int>(std::floor(radius / dr));
      if (bin >= radial_bins) {
        bin = radial_bins - 1;
      }
      if (bin < 0) {
        bin = 0;
      }
      rcount[static_cast<std::size_t>(bin)] += 1.0;
      rmass[static_cast<std::size_t>(bin)] += sample.mass;
    }
    const double theta = angleOf(y, x);
    const int abin = binOf((theta + kPi) / (2.0 * kPi), azimuthal_bins);
    tcount[static_cast<std::size_t>(abin)] += 1.0;
    tmass[static_cast<std::size_t>(abin)] += sample.mass;
  }
  finish(out.radial, rcount, rmass);
  finish(out.azimuthal, tcount, tmass);
  return out;
}

HostGuest hostGuest(const Cloud &cloud, int guest_type, int framework_type, double cutoff) {
  if (guest_type <= 0) {
    throw std::out_of_range("guest type must be positive");
  }
  if (framework_type == guest_type) {
    throw std::invalid_argument("framework type must differ from the guest type");
  }
  if (!std::isfinite(cutoff) || cutoff < 0.0) {
    throw std::out_of_range("cutoff must be non-negative");
  }
  const Cell cell = cellOf(cloud);
  requireCell(cell);
  std::vector<int> guests;
  std::vector<int> frames;
  for (int i = 0; i < static_cast<int>(cloud.pts.size()); ++i) {
    const int kind = cloud.pts[static_cast<std::size_t>(i)].type;
    if (kind == guest_type) {
      guests.push_back(i);
    } else if (framework_type == 0 || kind == framework_type) {
      frames.push_back(i);
    }
  }
  std::vector<int> types;
  for (int index : frames) {
    types.push_back(cloud.pts[static_cast<std::size_t>(index)].type);
  }
  std::sort(types.begin(), types.end());
  types.erase(std::unique(types.begin(), types.end()), types.end());
  HostGuest out;
  out.framework_types = types;
  const std::size_t nguest = guests.size();
  const std::size_t ntype = types.size();
  out.guest = guests;
  out.nearest.assign(nguest, -1);
  out.nearest_type.assign(nguest, 0);
  out.distance.assign(nguest, -1.0);
  out.contacts.assign(nguest * ntype, 0);
  const double cut2 = cutoff * cutoff;
  for (std::size_t g = 0; g < nguest; ++g) {
    const auto &guest = cloud.pts[static_cast<std::size_t>(guests[g])];
    double best = 0.0;
    bool have = false;
    for (int index : frames) {
      const auto &site = cloud.pts[static_cast<std::size_t>(index)];
      const auto delta =
          gen::triclinicMinImage(cloud, guest.x, guest.y, guest.z, site.x, site.y, site.z);
      const double r2 = delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2];
      if (!have || r2 < best) {
        have = true;
        best = r2;
        out.nearest[g] = index;
        out.nearest_type[g] = site.type;
        out.distance[g] = std::sqrt(r2);
      }
      if (r2 <= cut2) {
        const auto found = std::find(types.begin(), types.end(), site.type);
        if (found != types.end()) {
          const std::size_t column = static_cast<std::size_t>(found - types.begin());
          out.contacts[g * ntype + column] += 1;
        }
      }
    }
  }
  return out;
}

} // namespace pore
