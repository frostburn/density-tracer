# density-tracer
Ray tracer focused on self-illuminating density clouds and dispersive optics

## Compiling
```
mkdir build; cd build
cmake ..
make
```

## Testing
```
make test
make CTEST_OUTPUT_ON_FAILURE=1 test
```

## Running
```
./bin/render > /tmp/out.ppm
```

## Compiling a debug build
```
mkdir debug; cd debug
cmake -DCMAKE_BUILD_TYPE=Debug ..
make
```
