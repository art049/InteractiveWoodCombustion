CIBLE = main
SRCS =  Main.cpp Camera.cpp Mesh.cpp GLProgram.cpp GLShader.cpp GLError.cpp LightSource.cpp Tree.cpp
OPENGL_PATH = /usr/lib/nvidia-367 # change this for your own environment
LIBS = -L$(OPENGL_PATH) -lglut -lGLU -lGL -lGLEW -lm -lpthread

CC = g++
CPP = g++

FLAGS = -std=c++11 -O2 -Wall -g -pthread
LDFLAGS = -std=c++11
CFLAGS = $(FLAGS)
CXXFLAGS = $(FLAGS)

OBJS = $(SRCS:.cpp=.o)   

$(CIBLE): $(OBJS)
	g++ $(LDFLAGS) -o $(CIBLE) $(OBJS) $(LIBS)
clean:
	rm -f  *~  $(CIBLE) $(OBJS)

-include $(subst .cpp,.d,$(SRCS))
%.d : %.cpp
	$(CC) $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -MM -MF $@ -MP $<



