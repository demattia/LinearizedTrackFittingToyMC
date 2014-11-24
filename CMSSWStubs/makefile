# To get the ROOT includes do
# root-config --cflags
CFLAGS=-pthread -stdlib=libc++ -std=c++11 -m64 -I/Users/demattia/root/include
# For the libaries to link do
# root-config --libs
LDFLAGS=-L/Users/demattia/root/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lpthread -stdlib=libc++ -lm -ldl

all: buildMatrix

buildMatrix: main.o L1TrackTriggerTree.o MatrixBuilder.o
	clang $(LDFLAGS) main.o L1TrackTriggerTree.o MatrixBuilder.o -o buildMatrix

main.o: main.cpp
	clang $(CFLAGS) -c main.cpp

L1TrackTriggerTree.o: L1TrackTriggerTree.cpp
	clang $(CFLAGS) -c L1TrackTriggerTree.cpp

MatrixBuilder.o: MatrixBuilder.cpp
	clang $(CFLAGS) -c MatrixBuilder.cpp

clean:
	rm -rf *o main
