bbdirectory = ./harp

all: render_model_c #render_model_cuda

render_model_c: $(bbdirectory)/render_model.c $(bbdirectory)/render_model.h
	# mkdir -p lib
	gcc -shared -fPIC -O3 -o $(bbdirectory)/render_model.so $(bbdirectory)/render_model.c -I$(bbdirectory) -lm

# render_model_cuda: $(bbdirectory)/render_model_cuda.cu $(bbdirectory)/render_model_cuda.h
# 	nvcc --relocatable-device-code=true --compiler-options '-fPIC' -D_FORCE_INLINES -Xcompiler "-O3" -o $(bbdirectory)/render_model_cuda.o -c $(bbdirectory)/render_model_cuda.cu -I$(bbdirectory)
# 	nvcc --shared --relocatable-device-code=true --compiler-options '-fPIC' -D_FORCE_INLINES -Xcompiler "-O3" -o $(bbdirectory)/render_model_cuda.so $(bbdirectory)/render_model_cuda.o -I$(bbdirectory)
# 	rm $(bbdirectory)/render_model_cuda.o

clean:
	rm $(bbdirectory)/render_model.so
