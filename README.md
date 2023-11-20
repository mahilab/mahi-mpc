# Required Instructions

## IPOPT
Get ipopt from [here](https://github.com/coin-or/Ipopt/releases/tag/releases%2F3.14.3), and extract zip files. Place the whole folder into `C:/dev/coin-or`, so that it looks like
`C:/dev/coin-or/Ipopt-3.14.3-win64-msvs2019-md/`
Add the bin and lib directories from this directory to your system path.

## gcc
follow the installation instructions [here](https://www.msys2.org/) to get MSYS2. Going through step 7 will make gcc available. Note: It looks like instructions have changed since I first posted this. Don't click Run MSYS2 now in the installation helper. Instead, press the windows key and type in `MSYS2`. Make sure to open the one that has the name `MSYS2 MINGW64`. Then instead of the command it tells you to use to install gcc (`mingw-w64-ucrt-x86_64-gcc`), remove the `ucrt` and instead run 

```pacman -S mingw-w64-ucrt-x86_64-gcc```.

Add the directory `C:/msys64/mingw64/bin` to your path (check to make sure that folder has gcc in it).

## Summary
After this process, you should have 3 new directories on your path:
```
C:/dev/coin-or/Ipopt-3.14.3-win64-msvs2019-md/bin
C:/dev/coin-or/Ipopt-3.14.3-win64-msvs2019-md/lib
C:/msys64/mingw64/bin
```

If you need to know how to add something to your path in Windows, you can hit the windows key, search `path`, click on `Edit the system environment variables`, click on Environment Variables at the bottom, double click path (either under user variables or system variables). Then click `New` and enter the path you want to add.
