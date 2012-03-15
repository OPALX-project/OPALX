# serial version

#CC = mpiCC -DnoUSE_PBE -I$(IPL_ROOT)/src/pbe/include
CC = $(MYCXX) -Wall -DnoUSE_PBE -I$(IPL_ROOT)/src/pbe/include -I$(IPL_ROOT)/src/expde/2.1/program/source
CFLAGS = -Wno-deprecated  -ftemplate-depth-60 -g -funroll-loops -fstrict-aliasing  -DEXPDEPARALLEL -DCOMP_GNU -DWITH_BRICK
LDFLAGS = -static -lm -L../lib -lexpde

SOURCE = source
S_GRID = $(SOURCE)/grid
S_EXT  = $(SOURCE)/extemp
S_LIB  = $(SOURCE)/math_lib
S_DOM  = $(SOURCE)/domain
S_BAS  = $(SOURCE)/basic
S_FOR  = $(SOURCE)/formulas
S_EVP  = $(SOURCE)/evpar
S_IND  = $(SOURCE)/indices

H_EXTS  = $(S_EXT)/varcon.h \
          $(S_EXT)/variable.h

H_IND  =  $(S_IND)/index.h \

H_EXT  =  $(S_EXT)/variable.h \
          $(S_EXT)/variabl2.h \
          $(S_EXT)/opera.h \
          $(S_EXT)/operav.h \
          $(S_EXT)/op_conv.h \
          $(S_EXT)/mg_op.h \
          $(S_EXT)/res_op.h \
          $(S_EXT)/res_opv.h \
          $(S_EXT)/sim.h \
          $(S_EXT)/eleme.h \
          $(S_EXT)/extra.h \
          $(S_EXT)/array.h \
          $(S_EXT)/varcon.h \
	$(S_EXT)/partint.h \
          $(S_EVP)/evpar.h

H_DOE  =  $(S_DOM)/d_exam.h

H_GRID =  $(S_GRID)/gridbase.h  \
          $(S_GRID)/gpar.h \
          $(S_GRID)/input.h \
          $(S_GRID)/parallel.h \
          $(S_EVP)/evpar.h \
          $(S_GRID)/grid.h \
          $(S_GRID)/sto_man.h \
          $(S_GRID)/hash.h \
          $(S_GRID)/mgcoeff.h \
          $(S_GRID)/SSten.h \
          $(SOURCE)/access.h \
          $(S_GRID)/listint/listint1.h \
          $(S_GRID)/listint/listint2.h \
          $(S_GRID)/listint/listint3.h \
          $(S_GRID)/listint/listint4.h \
          $(S_GRID)/listnb/listnb1.h \
          $(S_GRID)/listnb/listnb2.h \
          $(S_GRID)/listnb/listnb3.h \
          $(S_GRID)/listnb/listnb4.h \
          $(S_GRID)/listbo/listbo1.h \
          $(S_GRID)/listpa/listpa1.h \
          $(S_GRID)/listce/listce1.h \
          $(S_GRID)/listcf/listcf1.h \
          $(S_GRID)/listcf/listcf2.h

H_FOR  =  $(S_FOR)/boundy.h
H_FO2  =  $(S_FOR)/diffop.h \
          $(S_FOR)/diffopc.h \
          $(S_FOR)/diffopcv.h \
          $(S_FOR)/diffopv.h
H_FO3  =  $(S_FOR)/subdiv.h \
          $(S_FOR)/subdivI.h \
          $(S_FOR)/subdivII.h \
          $(S_FOR)/subdiIII.h \
          $(S_FOR)/subdivIV.h
H_FO4  =  $(S_FOR)/loc_sten.h

H_DOM  =  $(S_DOM)/domain.h \
          $(S_DOM)/sca_dom.h
H_LIB  =  $(S_LIB)/math_lib.h

H_BAS  =  $(S_BAS)/basic.h \
          $(SOURCE)/access.h \
          $(SOURCE)/extpde.h \
          $(SOURCE)/paramete.h \
          $(SOURCE)/parser.h \
          $(SOURCE)/dummympi.h

OBJ = input.o \
      evpar.o \
      gpar.o \
      parallel.o \
      communi.o \
      commbos.o \
      mgcoeff.o \
      gridbase.o \
      labbo.o \
      hash.o \
      grid.o \
      ggen.o \
      printavs.o \
      printpa.o \
      p_opendx.o \
      pridxpa.o \
      domain.o \
      sca_dom.o \
      boundy.o \
      subdiv.o \
      loc_sten.o \
      res_sten.o \
      variable.o \
      index.o \
      neigh.o \
      varcon.o \
      eleme.o \
      opera.o \
      operav.o \
      op_conv.o \
      mg.o \
      opoth.o \
      math_lib.o \
	partint.o

all: $(OBJ)
	make lib
#make run
#	$(CC) $(CFLAGS) -o run $(OBJ) -lm 

run: main.o
	$(CC) $(CFLAGS) -o run main.o $(LDFLAGS)

partint.o : $(S_EXT)/partint.cc $(H_BAS) $(H_DOM) $(H_DOE) $(H_GRID) $(H_FOR) $(H_FO2) $(H_FO3) $(H_FO4) $(H_LIB) $(SOURCE)/abbrevi.h
	$(CC) -I$(EXPDE_ROOT)/program/source $(CFLAGS) -c $(S_EXT)/partint.cc

evpar.o : $(S_EVP)/evpar.cc $(H_BAS) $(H_DOM) $(H_DOE) $(H_GRID) $(H_FOR) $(H_FO2) $(H_FO3) $(H_FO4) $(H_LIB) $(SOURCE)/abbrevi.h
	$(CC) $(CFLAGS) -c $(S_EVP)/evpar.cc
parallel.o : $(S_GRID)/parallel.cc $(H_BAS) $(H_GRID) $(H_LIB) $(H_GRID) $(SOURCE)/abbrevi.h $(S_GRID)/parallel.h
	$(CC) $(CFLAGS) -c $(S_GRID)/parallel.cc
grid.o : $(S_GRID)/grid.cc $(H_BAS) $(H_DOM) $(H_DOE) $(H_GRID) $(H_FOR) $(H_FO2) $(H_FO3) $(H_FO4) $(H_LIB) $(SOURCE)/abbrevi.h
	$(CC) $(CFLAGS) -c $(S_GRID)/grid.cc
mgcoeff.o : $(S_GRID)/mgcoeff.cc $(H_BAS) $(H_DOM) $(H_DOE) $(H_GRID) $(H_FOR) $(H_FO2) $(H_FO3) $(H_FO4) $(H_LIB) $(SOURCE)/abbrevi.h
	$(CC) $(CFLAGS) -c $(S_GRID)/mgcoeff.cc
domain.o : $(S_DOM)/domain.cc $(H_BAS) $(H_DOM) $(H_GRID)  $(H_LIB) $(SOURCE)/abbrevi.h
	$(CC) $(CFLAGS) -c $(S_DOM)/domain.cc
sca_dom.o : $(S_DOM)/sca_dom.cc $(H_BAS) $(H_DOM) $(H_GRID)  $(H_LIB) $(SOURCE)/abbrevi.h
	$(CC) $(CFLAGS) -c $(S_DOM)/sca_dom.cc
boundy.o : $(S_FOR)/boundy.cc $(H_LIB) $(H_FOR) $(H_FO3) $(SOURCE)/abbrevi.h
	$(CC) $(CFLAGS) -c $(S_FOR)/boundy.cc
subdiv.o : $(S_FOR)/subdiv.cc $(H_LIB) $(H_BAS) $(H_FOR) $(H_FO3) $(SOURCE)/abbrevi.h
	$(CC) $(CFLAGS) -c $(S_FOR)/subdiv.cc
loc_sten.o : $(S_FOR)/loc_sten.cc $(H_LIB) $(H_FO4) $(SOURCE)/abbrevi.h
	$(CC) $(CFLAGS) -c $(S_FOR)/loc_sten.cc
res_sten.o : $(S_FOR)/res_sten.cc $(H_LIB) $(H_FO4) $(SOURCE)/abbrevi.h
	$(CC) $(CFLAGS) -c $(S_FOR)/res_sten.cc
math_lib.o : $(S_LIB)/math_lib.cc $(H_LIB)
	$(CC) $(CFLAGS) -c $(S_LIB)/math_lib.cc
variable.o : $(S_EXT)/variable.cc $(H_DOM)  $(H_EXTS) $(H_BAS) $(H_LIB) $(H_FOR) $(H_FO2)  $(H_FO3) $(H_GRID) $(SOURCE)/abbrevi.h
	$(CC) $(CFLAGS) -c $(S_EXT)/variable.cc
index.o : $(S_IND)/index.cc $(H_DOM)  $(H_EXTS) $(H_BAS) $(H_LIB) $(H_FOR) $(H_IND) $(H_FO2) $(H_FO3) $(H_EXT) $(H_GRID) $(SOURCE)/abbrevi.h
	$(CC) $(CFLAGS) -c $(S_IND)/index.cc
neigh.o : $(S_IND)/neigh.cc $(H_DOM)  $(H_EXTS) $(H_BAS) $(H_LIB) $(H_FOR) $(H_IND) $(H_FO2) $(H_FO3) $(H_EXT) $(H_GRID) $(SOURCE)/abbrevi.h
	$(CC) $(CFLAGS) -c $(S_IND)/neigh.cc
varcon.o : $(S_EXT)/varcon.cc $(H_DOM)  $(H_EXTS) $(H_BAS) $(H_LIB) $(H_FOR) $(H_FO2)  $(H_FO3) $(H_GRID) $(SOURCE)/abbrevi.h
	$(CC) $(CFLAGS) -c $(S_EXT)/varcon.cc
eleme.o : $(S_EXT)/eleme.cc $(H_DOM)  $(H_EXT) $(H_BAS) $(H_LIB) $(H_FOR) $(H_FO2)  $(H_FO3) $(H_GRID) $(SOURCE)/abbrevi.h
	$(CC) $(CFLAGS) -c $(S_EXT)/eleme.cc
opera.o : $(S_EXT)/opera.cc $(H_DOM)  $(H_EXT) $(H_BAS) $(H_LIB) $(H_FOR) $(H_FO2)  $(H_FO3) $(H_GRID) $(SOURCE)/abbrevi.h
	$(CC) $(CFLAGS) -c $(S_EXT)/opera.cc
opoth.o : $(S_EXT)/opoth.cc $(H_DOM)  $(H_EXT) $(H_BAS) $(H_LIB) $(H_FOR) $(H_FO2)  $(H_FO3) $(H_GRID) $(SOURCE)/abbrevi.h
	$(CC) $(CFLAGS) -c $(S_EXT)/opoth.cc
operav.o : $(S_EXT)/operav.cc $(H_DOM)  $(H_EXT) $(H_BAS) $(H_LIB) $(H_FOR) $(H_FO2)  $(H_FO3) $(H_GRID) $(SOURCE)/abbrevi.h
	$(CC) $(CFLAGS) -c $(S_EXT)/operav.cc
op_conv.o : $(S_EXT)/op_conv.cc $(H_DOM)  $(H_EXT) $(H_BAS) $(H_LIB) $(H_FOR) $(H_FO2)  $(H_FO3) $(H_GRID) $(SOURCE)/abbrevi.h
	$(CC) $(CFLAGS) -c $(S_EXT)/op_conv.cc
mg.o : $(S_EXT)/mg.cc $(H_DOM)  $(H_EXT) $(H_BAS) $(H_LIB) $(H_FOR) $(H_FO2)  $(H_FO3) $(H_GRID) $(SOURCE)/abbrevi.h
	$(CC) $(CFLAGS) -c $(S_EXT)/mg.cc
gridbase.o : $(S_GRID)/gridbase.cc $(H_DOM) $(H_BAS) $(H_FOR) $(H_FO3) $(H_GRID)  $(H_LIB) $(SOURCE)/abbrevi.h
	$(CC) $(CFLAGS) -c $(S_GRID)/gridbase.cc
gpar.o : $(S_GRID)/gpar.cc $(H_DOM) $(H_BAS) $(H_FOR) $(H_FO3) $(H_GRID)  $(H_LIB) $(SOURCE)/abbrevi.h
	$(CC) $(CFLAGS) -c $(S_GRID)/gpar.cc
communi.o : $(S_EVP)/communi.cc $(H_DOM) $(H_BAS) $(H_FOR) $(H_FO3) $(H_GRID)  $(H_LIB) $(SOURCE)/abbrevi.h
	$(CC) $(CFLAGS) -c $(S_EVP)/communi.cc
commbos.o : $(S_EVP)/commbos.cc $(H_DOM) $(H_BAS) $(H_FOR) $(H_FO3) $(H_GRID)  $(H_LIB) $(SOURCE)/abbrevi.h
	$(CC) $(CFLAGS) -c $(S_EVP)/commbos.cc
labbo.o : $(S_GRID)/labbo.cc $(H_DOM) $(H_BAS) $(H_FOR) $(H_FO3) $(H_GRID)  $(H_LIB) $(SOURCE)/abbrevi.h
	$(CC) $(CFLAGS) -c $(S_GRID)/labbo.cc
hash.o : $(S_GRID)/hash.cc $(H_DOM) $(H_BAS) $(H_FOR) $(H_FO3) $(H_GRID)  $(H_LIB) $(SOURCE)/abbrevi.h
	$(CC) $(CFLAGS) -c $(S_GRID)/hash.cc
ggen.o : $(S_GRID)/ggen.cc $(H_DOM) $(H_BAS) $(H_FOR) $(H_FO3) $(H_GRID)  $(H_LIB) $(SOURCE)/abbrevi.h
	$(CC) $(CFLAGS) -c $(S_GRID)/ggen.cc
printavs.o : $(S_GRID)/printavs.cc $(H_BAS) $(H_GRID) $(H_FOR)  $(H_FO3) $(H_LIB) $(SOURCE)/abbrevi.h
	$(CC) $(CFLAGS) -c $(S_GRID)/printavs.cc
input.o : $(S_GRID)/input.cc $(H_BAS) $(H_GRID) $(H_FOR)  $(H_FO3) $(H_LIB) $(SOURCE)/abbrevi.h
	$(CC) $(CFLAGS) -c $(S_GRID)/input.cc
printpa.o : $(S_GRID)/printpa.cc $(H_BAS) $(H_GRID) $(H_FOR)  $(H_FO3) $(H_LIB) $(SOURCE)/abbrevi.h
	$(CC) $(CFLAGS) -c $(S_GRID)/printpa.cc
p_opendx.o : $(S_GRID)/p_opendx.cc $(H_BAS) $(H_GRID) $(H_FOR)  $(H_FO3) $(H_LIB) $(SOURCE)/abbrevi.h
	$(CC) $(CFLAGS) -c $(S_GRID)/p_opendx.cc
pridxpa.o : $(S_GRID)/pridxpa.cc $(H_BAS) $(H_GRID) $(H_FOR)  $(H_FO3) $(H_LIB) $(SOURCE)/abbrevi.h
	$(CC) $(CFLAGS) -c $(S_GRID)/pridxpa.cc


# rm main.o because we do not want a main in the library

lib: 	$(OBJ)
	mkdir -p ../lib
	ar rv ../lib/libexpde.a *.o

clean :
	rm run $(OBJ)



