import os
import subprocess

# 定义比率列表
ratios = [1,0.1,0.01,0.001]

# 遍历每个比率
for ratio in ratios:
    # 读取parameter.py的内容
    with open("disk_properties.py", "r") as f:
        lines = f.readlines()

    # 替换包含"ratio ="的行
    for i, line in enumerate(lines):
        if "ratio =" in line:
            lines[i] = f"ratio = {ratio} \n"

    # 将新的内容写入到parameter.py
    with open("disk_properties.py", "w") as f:
        f.writelines(lines)

    import pdb;pdb.set_trace()
    # 运行runcpd.py
    subprocess.run(["python", "runcpd.py"])
