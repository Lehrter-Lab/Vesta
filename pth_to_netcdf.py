import numpy as np
from netCDF4 import Dataset

def init_netcdf(outfile):
    root = Dataset(outfile, "w", format="NETCDF4")
    root.createDimension("point", None)  # unlimited
    comp = dict(zlib=True, complevel=4, shuffle=True)
    
    root.createVariable("pid", "i4", ("point",), zlib=False)
    for name in ["lon", "lat", "depth", "time"]:
        root.createVariable(name, "f4", ("point",), **comp)

    return root

def append_chunk(root, buffer, start_index):
    n     = len(buffer)
    pid   = np.array([b[0] for b in buffer], dtype=np.int32)
    x     = np.array([b[1] for b in buffer], dtype=np.float32)
    y     = np.array([b[2] for b in buffer], dtype=np.float32)
    depth = np.array([b[3] for b in buffer], dtype=np.float32)
    t     = np.array([b[4] for b in buffer], dtype=np.float64)

    end_index = start_index + n
    root["pid"][start_index:end_index]   = pid
    root["lon"][start_index:end_index]   = x
    root["lat"][start_index:end_index]   = y
    root["depth"][start_index:end_index] = depth
    root["time"][start_index:end_index]  = t

    return end_index

def pth_to_netcdf(pth_file, output_nc, chunk_size=100000):
    root = init_netcdf(output_nc)
    buffer = []
    i = 0

    with open(pth_file) as f:
        while True:
            header = f.readline()
            if not header:
                break
            parts = header.strip().split()
            if len(parts) != 2:
                continue
            try:
                time_sec, n_particles = float(parts[0]), int(parts[1])
            except ValueError:
                continue

            for _ in range(n_particles):
                line = f.readline()
                if not line:
                    raise EOFError("Unexpected EOF while reading particle data")
                vals = line.strip().split()
                if len(vals) < 4:
                    raise ValueError(f"Malformed particle line: {line}")

                buffer.append((int(vals[0]), float(vals[1]), float(vals[2]), float(vals[3]), time_sec))

                if len(buffer) >= chunk_size:
                    i = append_chunk(root, buffer, i)
                    buffer = []

        if buffer:
            i = append_chunk(root, buffer, i)

    root.close()
    print(f"Finished writing {i} particles to {output_nc}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 3:
        print("Usage: python pth_to_netcdf.py particle.pth output.nc [chunk_size]")
        sys.exit(1)
    pth_file, output_nc = sys.argv[1], sys.argv[2]
    #pth_file, output_nc = 'particle.pth', 'output.nc'
    chunk_size = int(sys.argv[3]) if len(sys.argv) == 4 else 100000
    pth_to_netcdf(pth_file, output_nc, chunk_size)

