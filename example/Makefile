include ../paths.mk
CPP_FILES := $(wildcard *.cpp)
TARGETS := $(addprefix bin/,$(notdir $(CPP_FILES:.cpp=)))
LD_FLAGS := -L../bin -lcollisions
CC_FLAGS := -g -Wall -O2

all: $(TARGETS)

bin/%: obj/%.o
	g++ -o $@ $^ $(LD_FLAGS)

obj/%.o: %.cpp
	g++ $(CC_FLAGS) -I$(EIGEN_DIR) -I../include -c -o $@ $<

clean:
	rm $(TARGETS)
