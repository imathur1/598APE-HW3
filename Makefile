FUNC := clang++
copt := -c 
OBJ_DIR := ./bin/
FLAGS := -O3 -lm -g -Werror -lprofiler

# Enable profiling via a macro
ifndef PROFILE
	PROFILE := 0
endif
FLAGS += -DPROFILE=$(PROFILE)

CPP_FILES := $(wildcard src/*.cpp)
OBJ_FILES := $(addprefix $(OBJ_DIR),$(notdir $(CPP_FILES:.cpp=.obj)))

all:
	$(FUNC) ./main.cpp -o ./main.exe $(FLAGS)

clean:
	rm -f ./*.exe

view-profile:
	~/go/bin/pprof -http "0.0.0.0:8080" ./main.exe ./my_profile.prof

.PHONY: all clean view-profile