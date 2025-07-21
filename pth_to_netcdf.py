import numpy as np
import xarray as xr

def count_total_particles(filepath):
    total = 0
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) != 2:
                continue
            try:
                _ = float(parts[0])
                n = int(parts[1])
            except ValueError:
                continue
            total += n
    return total

def read_pth_with_depth(filepath):
    total_points = count_total_particles(filepath)
    x     = np.empty(total_points, dtype=np.float64)
    y     = np.empty(total_points, dtype=np.float64)
    t     = np.empty(total_points, dtype=np.float64)
    pid   = np.empty(total_points, dtype=np.int32)
    depth = np.empty(total_points, dtype=np.float32)
    
    with open(filepath) as f:
        i = 0
        while True:
            header = f.readline()
            if not header:
                break
            parts = header.strip().split()
            if len(parts) != 2:
                continue
            try:
                time_sec    = float(parts[0])
                n_particles = int(parts[1])
            except ValueError:
                continue
            
            for _ in range(n_particles):
                line = f.readline()
                if not line:
                    raise EOFError("Unexpected EOF while reading particle data")
                vals = line.strip().split()
                if len(vals) < 4:
                    raise ValueError(f"Malformed particle line: {line}")
                pid[i]   = int(vals[0])
                x[i]     = float(vals[1])
                y[i]     = float(vals[2])
                depth[i] = float(vals[3])
                t[i]     = time_sec
                i += 1
    
    if i != total_points:
        raise RuntimeError(f"Read {i} particles but expected {total_points}")
    return pid, t, x, y, depth

def write_to_netcdf(outfile, pid, t, x, y, depth):
    ds = xr.Dataset(
        {
            "time": (["point"], t),
            "pid": (["point"], pid),
            "longitude": (["point"], x),
            "latitude": (["point"], y),
            "depth": (["point"], depth),
        },
        coords={"point": np.arange(len(pid))},
        attrs={"description": "Raw particle tracking data from .pth file"},
    )
    comp     = dict(zlib=True, complevel=4)
    encoding = {var: comp for var in ds.data_vars}
    ds.to_netcdf(outfile, encoding=encoding)
    print(f"Saved raw particle data to {outfile}")

def main(pth_file, output_nc):
    pid, t, x, y, depth = read_pth_with_depth(pth_file)
    write_to_netcdf(output_nc, pid, t, x, y, depth)

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python pth_to_netcdf.py particle.pth output.nc")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
