all:	library programs util gammel subopt perl

library:
	cd lib;	$(MAKE)

programs:
	cd Progs; $(MAKE)

util:
	cd Utils; $(MAKE)

gammel:
	cd Cluster; $(MAKE)

subopt:
	cd Subopt; $(MAKE)

perl:
	cd Perl; $(MAKE) -f Makefile.old Makefile; $(MAKE)

install:	all
	cd Progs; $(MAKE) install
	cd Cluster; $(MAKE) install
	cd Subopt; $(MAKE) install
	cd Utils; $(MAKE) install 
	echo "type  cd Perl; $(MAKE) install   to install the Perl modules" 

clean:	
	cd lib; $(MAKE) clean
	cd Progs; $(MAKE) clean
	cd Cluster; $(MAKE) clean
	cd Subopt; $(MAKE) clean
	cd Utils; $(MAKE) clean
	cd Perl; $(MAKE) clean
