# Makefile for Fortran + MATLAB Engine
FC = gfortran
CC = gcc

# MATLAB路径（需要根据你的安装路径调整）
MATLAB_ROOT = /Applications/MATLAB_R2024b.app
MATLAB_INC = $(MATLAB_ROOT)/extern/include
MATLAB_LIB = $(MATLAB_ROOT)/bin/maci64

# 编译选项
FFLAGS = -O2 -Wall
CFLAGS = -O2 -Wall
LDFLAGS = -L$(MATLAB_LIB) -leng -lmx -lm

# 目标文件
TARGET = main_with_matlab
SOURCE = main_with_matlab.f90

# 默认目标
all: $(TARGET)

# 编译规则
$(TARGET): $(SOURCE)
	$(FC) $(FFLAGS) -I$(MATLAB_INC) -o $@ $< $(LDFLAGS)

# 清理
clean:
	rm -f $(TARGET) *.o *.mod

# 运行
run: $(TARGET)
	./$(TARGET)

# 帮助
help:
	@echo "可用的目标："
	@echo "  all     - 编译程序"
	@echo "  clean   - 清理编译文件"
	@echo "  run     - 编译并运行程序"
	@echo "  help    - 显示此帮助信息"

.PHONY: all clean run help 