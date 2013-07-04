#!/bin/bash
#

WIN_INSTALLER_LOG=win_installer.log

export ac_cv_func_realloc_0_nonnull=yes
export ac_cv_func_malloc_0_nonnull=yes



case "$1" in

arch-mingw64) echo "Making Windows Installer using Arch Linux mingw32 installation of mingw-w64"
              echo -ne "...making 32bit version..."
              ./configure --host=i486-mingw32 --without-perl --without-forester --with-cluster &>> ${WIN_INSTALLER_LOG}
              make clean &>> ${WIN_INSTALLER_LOG}
              make &>> ${WIN_INSTALLER_LOG}
              makensis win_installer_archlinux_i486.nsi &>> ${WIN_INSTALLER_LOG}
              echo -ne " done\n"
              ;;

arch)         echo "Making Windows Installer using Arch Linux mingw32 installation of mingw"
              echo -ne "...making 32bit version..."
              ./configure --host=i486-mingw32 --without-perl --without-forester --with-cluster &>> ${WIN_INSTALLER_LOG}
              make clean &>> ${WIN_INSTALLER_LOG}
              make &>> ${WIN_INSTALLER_LOG}
              makensis win_installer_archlinux_i486.nsi &>> ${WIN_INSTALLER_LOG}
              echo -ne " done\n"
              ;;

*)            echo "Making Windows Installer using Fedora installation of mingw64"
              echo -ne "...making 32bit version..."
              mingw32-configure --without-perl --without-forester --with-cluster &>> ${WIN_INSTALLER_LOG}
              make clean  &>> ${WIN_INSTALLER_LOG}
              make &>> ${WIN_INSTALLER_LOG}
              makensis win_installer_fedora_i686.nsi >> ${WIN_INSTALLER_LOG}
              echo -ne " done\n"
              echo -ne "...making 64bit version..."
              mingw64-configure --without-perl --without-forester --with-cluster &>> ${WIN_INSTALLER_LOG}
              make clean  &>> ${WIN_INSTALLER_LOG}
              make &>> ${WIN_INSTALLER_LOG}
              makensis win_installer_fedora_x86_64.nsi &>> ${WIN_INSTALLER_LOG}
              echo -ne "done\n"
              ;;

esac


