# Makefile for k_solver Fortran program
FC = gfortran
FFLAGS = -O2 -Wall -Wextra -std=f2008
TARGET = k_solver
SOURCES = k_solver.f90
OBJECTS = $(SOURCES:.f90=.o)

# 默认目标
all: $(TARGET)

# 编译主程序
$(TARGET): $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $^

# 编译目标文件
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

# 清理
clean:
	rm -f $(OBJECTS) $(TARGET)

# 运行程序
run: $(TARGET)
	./$(TARGET)

# 调试版本
debug: FFLAGS += -g -fcheck=all
debug: $(TARGET)

# 显示帮助
help:
	@echo "可用的目标:"
	@echo "  all     - 编译程序 (默认)"
	@echo "  clean   - 清理编译文件"
	@echo "  run     - 编译并运行程序"
	@echo "  debug   - 编译调试版本"
	@echo "  help    - 显示此帮助信息"

.PHONY: all clean run debug help 