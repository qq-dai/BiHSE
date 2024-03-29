# Compile
Using cmake files to compile the source files. For instance,
```
cd biclique
cmake .
make
```

# Usage
To execute the code, you need to run the following executable files, which accept the following optional parameters:

- "-f": the running bipartite graph.

- "-l": The left size constraint.

- "-r": The right size constraint.

- "-k": the value of k for k-plex. 

- "-d": Selceted from 'two' or 'core'.

An running example:

```
./biclique/bin/MBC -f datas/moreno.txt -d core
```

In 'BCE_ploy.cpp', we also implement a polynomial delay algorithm for maximal biclique enumerations. To compile this file, just change the file name of 'BCE_poly.cpp' to 'BCE.cpp'.
