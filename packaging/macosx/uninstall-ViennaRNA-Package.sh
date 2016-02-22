#!/bin/bash -e
#
# Uninstall the ViennaRNA Package for MacOSX.
#


if [ ! -x "/usr/local/bin/RNAfold" ] ; then
  echo "The ViennaRNA Package doesn't appear to be installed. Aborting"
  exit 1
fi

echo "This will uninstall the ViennaRNA Package"
printf "Type 'yes' if you are sure you wish to continue: "
read response
if [ "$response" == "yes" ]
then

  # remove all files installed by this package
  for pkg in `pkgutil --packages | grep org.TBI.ViennaRNA`
  do
    # note, pkgutil lists files without root path, i.e.
    # usr/local/bin/RNAfold instead of /usr/local/bin/RNAfold
    # Thus, we need to add the root '/' ourselves here
    pkgutil --only-files --files $pkg | xargs -I {} sudo rm -f /{}
  done

  # forget that ViennaRNA Package was installed
  echo "Removing ViennaRNA Package from pkgutil database"
  pkgutil --packages | grep org.TBI.ViennaRNA | xargs -I {} sudo pkgutil --forget {}

  # cleanup some potential remnants
  sudo rm -rf /usr/local/share/ViennaRNA
  sudo rm -rf /usr/local/share/doc/ViennaRNA

  echo "Successfully uninstalled the ViennaRNA Package"
else
  echo "Aborted"
  exit 1
fi

exit 0
