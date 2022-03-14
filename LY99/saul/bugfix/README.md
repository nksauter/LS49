This script fails just before completion.

[reproducer.sh](./reproducer.sh).  The *.err file shows this: 
```terminate called after throwing an instance of 'std::runtime_error'
  what():  Kokkos allocation "device_Fhkl" is being deallocated after Kokkos::finalize was called

Traceback functionality not available``` 

To run the reproducer, follow this:
```
$ cd $MODULES
$ git clone https://gitlab.com/cctbx/psii_spread.git
$ mkdir ~/.cctbx.xfel
$ ln -s $PWD/psii_spread/merging ~/.cctbx.xfel/
```

