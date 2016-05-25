
INCLUDES =   -I./includes
CXXFLAGS =	-O2 -ggdb -Wall -fmessage-length=0 $(INCLUDES)



OBJS =		./src/gsm_decoder.o  ./src/extra_functions.o  ./src/sch.o  ./src/receiver_config.o

LIBS = -lgnuradio-blocks -lgnuradio-runtime -lboost_filesystem -lliquid

TARGET =	gsm_decoder

DEFAULT_INCLUDES = -I./includes

#SCH= sch.o

#$(SCH) : gcc -I./includes -c ./src/sch.c -o ./src/sch.o

#OBJ2 = ./src/sch.o

$(TARGET):	$(OBJS)
	$(CXX)  $(INCLUDES)  -o $(TARGET) $(OBJS)   $(LIBS) $(DEFAULT_INCLUDES)

all:	$(TARGET) 

clean:
	rm -f $(OBJS) $(TARGET)
