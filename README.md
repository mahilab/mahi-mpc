# Required Instructions

## IPOPT
Get ipopt from [here](https://github.com/coin-or/Ipopt/releases/tag/releases%2F3.14.3), and extract zip files. Place the whole folder into `C:/dev/coin-or`, so that it looks like
`C:/dev/coin-or/Ipopt-3.14.3-win64-msvs2019-md/`
Add the bin and lib directories from this directory to your system path.

## gcc
follow the installation instructions [here](https://www.msys2.org/) to get MSYS2. Going through step 7 will make gcc available. 
Add the directory `C:/msys64/mingw64/bin` to your path.

## Summary
After this process, you should have 3 new directories on your path:
```
C:/dev/coin-or/Ipopt-3.14.3-win64-msvs2019-md/bin
C:/dev/coin-or/Ipopt-3.14.3-win64-msvs2019-md/lib
C:/msys64/mingw64/bin
```

If you need to know how to add something to your path in Windows, you can hit the windows key, search `path`, click on `Edit the system environment variables`, click on Environment Variables at the bottom, double click path (either under user variables or system variables). Then click `New` and enter the path you want to add.