
This repo contains code and plots for the paper Efficiently Computing Maximum Flows in Scale-Free Networks.

## Dependencies

- cmake
- C++17
- openMP
- boost
- googletest (a copy or symlink is expected in `code/tests/googletest`)
- [girgs](https://github.com/chistopher/girgs)

## Build

The code can be built using cmake. In the `code` directory do the following:
```
mkdir build
cd buid
cmake ..
make -j 8
```

## Data

The network data and raw results of experiements can be found in the release and should be unpacked to `code/data`.

