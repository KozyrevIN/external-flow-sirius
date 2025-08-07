# external-flow-sirius
Repo for the project on the flow around the body problem and calculating surface derivatives.
for run common mode:
cmake -S . -B build && cmake --build build --target main && ./build/main
for run debug mode with sanitizers:
cmake -S . -B build && cmake --build build --target debug && ./build/debug