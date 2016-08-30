#!/bin/bash
#

WIN_INSTALLER_LOG=win_installer.log

export ac_cv_func_realloc_0_nonnull=yes
export ac_cv_func_malloc_0_nonnull=yes

CONFIGURE_OPTIONS=" --without-swig \
                    --without-doc \
                    --without-forester \
                    --with-cluster \
                    --with-kinwalker \
                    --disable-tty-colors"

case "$1" in

arch64) echo "Making Windows Installer using Arch Linux mingw-w64 installation"
              echo -ne "...making 32bit version..."
              cd ..
              ./configure --host=i686-w64-mingw32 ${CONFIGURE_OPTIONS} >> packaging/${WIN_INSTALLER_LOG} 2>&1
              make clean >> packaging/${WIN_INSTALLER_LOG} 2>&1
              make >>${WIN_INSTALLER_LOG} 2>&1
              cd packaging
              makensis win_installer_archlinux_i686.nsi >>${WIN_INSTALLER_LOG} 2>&1
              echo -ne " done\n"
              echo -ne "...making 64bit version..."
              cd ..
              ./configure --host=x86_64-w64-mingw32 ${CONFIGURE_OPTIONS} >> packaging/${WIN_INSTALLER_LOG} 2>&1
              make clean >> packaging/${WIN_INSTALLER_LOG} 2>&1
              make >>${WIN_INSTALLER_LOG} 2>&1
              cd packaging
              makensis win_installer_archlinux_x86_64.nsi >>${WIN_INSTALLER_LOG} 2>&1
              echo -ne " done\n"
              ;;

arch)         echo "Making Windows Installer using Arch Linux mingw32 installation of mingw"
              echo -ne "...making 32bit version..."
              ./configure --host=i486-mingw32 ${CONFIGURE_OPTIONS} >>${WIN_INSTALLER_LOG} 2>&1
              make clean >>${WIN_INSTALLER_LOG} 2>&1
              make >>${WIN_INSTALLER_LOG} 2>&1
              makensis packaging/win_installer_archlinux_i486.nsi >>${WIN_INSTALLER_LOG} 2>&1
              echo -ne " done\n"
              ;;

*)            echo "Making Windows Installer using Fedora installation of mingw64"
              echo -ne "...making 32bit version..."
              mingw32-configure --without-perl ${CONFIGURE_OPTIONS} >>${WIN_INSTALLER_LOG} 2>&1
              make clean  >>${WIN_INSTALLER_LOG} 2>&1
              make >>${WIN_INSTALLER_LOG} 2>&1
              makensis packaging/win_installer_fedora_i686.nsi >> ${WIN_INSTALLER_LOG}
              echo -ne " done\n"
              echo -ne "...making 64bit version..."
              mingw64-configure --without-perl ${CONFIGURE_OPTIONS} >>${WIN_INSTALLER_LOG} 2>&1
              make clean  >>${WIN_INSTALLER_LOG} 2>&1
              make >>${WIN_INSTALLER_LOG} 2>&1
              makensis packaging/win_installer_fedora_x86_64.nsi >>${WIN_INSTALLER_LOG} 2>&1
              echo -ne "done\n"
              ;;

esac
