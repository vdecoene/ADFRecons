CPP := g++
CPPFLAGS := -O2 -Wno-write-strings
FC := gfortran
FFLAGS := -O2 -fno-second-underscore -w
LIBS := -L/usr/lib -lc -lm -lgfortran

all: lib/libport_i.a

clean:
	@rm -rf objs lib example

example: objs/example.o objs/port_i.o objs/port3f.o
	@$(CPP) $(CPPFLAGS) -o example objs/example.o objs/port_i.o objs/port3f.o $(LIBS)

objs/example.o: src/example.cxx objs
	@$(CPP) $(CPPFLAGS) -o objs/example.o -c src/example.cxx

objs/port_i.o: src/port_i.h src/port_i.cxx objs
	@$(CPP) $(CPPFLAGS) -o objs/port_i.o -c src/port_i.cxx

objs/port3f.o: port/port3f.f objs
	@rm -f port/port3f.f
	@cat port/*.f > port/port3f.f
	@$(FC) $(FFLAGS) -o objs/port3f.o -c port/port3f.f

lib/libport_i.a: objs/port_i.o objs/port3f.o
	@mkdir -p lib
	@ar rv lib/libport_i.a objs/port_i.o objs/port3f.o
	@ranlib lib/libport_i.a

objs:
	@mkdir -p objs
