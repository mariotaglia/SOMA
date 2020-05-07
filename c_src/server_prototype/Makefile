# Makefile for gtest examples

GOOGLE_TEST_LIB = gtest
GOOGLE_TEST_INCLUDE = /usr/local/include

G++ = mpicxx
G++_FLAGS = -c -Werror -I $(GOOGLE_TEST_INCLUDE)
LD_FLAGS = -L /usr/local/lib -l $(GOOGLE_TEST_LIB) -l pthread

OBJECTS = unittests.o server.o serv_util.o
TARGET = runtests

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(G++) -o $(TARGET) $(OBJECTS) $(LD_FLAGS)

unittests.o : unittests.cpp
	$(G++) $(G++_FLAGS) $<

server.o : server.c
	$(G++) $(G++_FLAGS) $<

serv_util.o : serv_util.c
	$(G++) $(G++_FLAGS) $<

clean:
	rm -f $(TARGET) $(OBJECTS)
                    
.PHONY: all clean
