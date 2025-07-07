import numpy as np
from scipy import special

# 测试相同的输入值
test_cases = [1+1j, 1000+1000j, 0.5+0.5j]

print("Python scipy.special.wofz 结果:")
print("=" * 50)

for test_z in test_cases:
    result = special.wofz(test_z)
    print(f"输入: {test_z}")
    print(f"wofz({test_z}) = {result}")
    print(f"实部: {result.real:.6e}")
    print(f"虚部: {result.imag:.6e}")
    print("-" * 30)

# 特别检查 wofz(1+1j) 的精确值
z = 1 + 1j
result = special.wofz(z)
print(f"\nwofz(1+1j) 的精确值:")
print(f"result = {result}")
print(f"实部: {result.real:.15e}")
print(f"虚部: {result.imag:.15e}") 