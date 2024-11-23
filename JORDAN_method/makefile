CXXFLAGS = -O3 -mfpmath=sse -fstack-protector-all -g \
           -W -Wall -Wextra -Wunused -Wcast-align \
           -Werror -pedantic -pedantic-errors -Wfloat-equal \
           -Wpointer-arith -Wformat-security -Wmissing-format-attribute \
           -Wformat=1 -Wwrite-strings -Wno-long-long \
           -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual \
           -Wno-suggest-attribute=format

TARGET = a.out

SRCS = main.cpp functions.h


all:
	g++ $(CXXFLAGS) $(SRCS) -o $(TARGET)
    # g++ -g $(SRCS) -o $(TARGET)
