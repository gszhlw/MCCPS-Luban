import random
import math
import copy

# 设置盒子中的原子个数
n_atoms = 25

# 设置执行蒙特卡洛移动次数
num_moves = 5000

# 设置模拟盒尺寸（单位：埃）
box_size = [ 15.0, 15.0, 15.0 ]

# 原子位移的最大值
max_translate = 0.5    # 单位：埃

# 模拟温度
temperature = 298.15   # k温度：开尔文

# 给出原子的Lennard Jones参数
# 用于氪原子的OPLS参数
sigma = 3.624         # 单位：埃
epsilon = 0.317       # 单位：kcal mol-1

# 创建数组以保存原子坐标
coords = []

# 盒中原子坐标随机生成
for i in range(0,n_atoms):

    coords.append( [random.uniform(0,box_size[0]), \
                    random.uniform(0,box_size[1]), \
                    random.uniform(0,box_size[2]) ] )


def make_periodic(x, box):
    """周期性边界条件（PBC）"""
    while (x < -0.5*box):
        x += box

    while (x > 0.5*box):
        x -= box

    return x


def wrap_into_box(x, box):
    """将盒中原子坐标包装进盒子（保证随机位移后不会移出盒子）"""
    while (x > box):
        x -= box

    while (x < 0):
        x += box

    return x


def print_pdb(move):
    """P对于特定移动打印PDB文件（可视化）"""
    filename = "output%000006d.pdb" % move

    FILE = open(filename, "w")

    FILE.write("CRYST1 %8.3f %8.3f %8.3f  90.00  90.00  90.00\n" % \
               (box_size[0], box_size[1], box_size[2]))

    for i in range(0,n_atoms):
        FILE.write("ATOM  %5d  Kr   Kr     1    %8.3f%8.3f%8.3f  1.00  0.00          Kr\n" % \
                   (i+1, coords[i][0], coords[i][1], coords[i][2]))
        FILE.write("TER\n")

    FILE.close()


# 计算原子能量
def calculate_energy():
    """计算传递原子的嫩能量（假设所有的原子对应的LJ的sigma值、epsilon值相同）"""

    # 循环所有原子对并计算其LJ能量

    total_energy = 0

    for i in range(0,n_atoms-1):
        for j in range(i+1, n_atoms):
            delta_x = coords[j][0] - coords[i][0]
            delta_y = coords[j][1] - coords[i][1]
            delta_z = coords[j][2] - coords[i][2]

            # 应用周期性边界条件
            delta_x = make_periodic(delta_x, box_size[0])
            delta_y = make_periodic(delta_y, box_size[1])
            delta_z = make_periodic(delta_z, box_size[2])

            # 计算原子间距离
            r = math.sqrt( (delta_x*delta_x) + (delta_y*delta_y) +
                           (delta_z*delta_z) )

            # E_LJ = 4*epsilon[ (sigma/r)^12 - (sigma/r)^6 ]
            e_lj = 4.0 * epsilon * ( (sigma/r)**12 - (sigma/r)**6 )

            total_energy += e_lj

    # 返回原子的总能量
    return total_energy


# 计算 kT
k_boltz = 1.987206504191549E-003  # kcal mol-1 K-1

kT = k_boltz * temperature

# 接受移动总数
naccept = 0

# 拒绝移动总数
nreject = 0

# 打印初始PDB文件
print_pdb(0)

# 执行所有移动
for move in range(1,num_moves+1):

    # 计算旧能量
    old_energy = calculate_energy()

    # 随机选择一个原子，即随机选择x和y之间的随机整数（包括x和y）
    atom = random.randint(0, n_atoms-1)

    # 保存旧坐标
    old_coords = copy.deepcopy(coords)

    # 执行移动——每个维度以delta值移动
    delta_x = random.uniform( -max_translate, max_translate )
    delta_y = random.uniform( -max_translate, max_translate )
    delta_z = random.uniform( -max_translate, max_translate )

    coords[atom][0] += delta_x
    coords[atom][1] += delta_y
    coords[atom][2] += delta_z

    # 将坐标包装回盒中
    coords[atom][0] = wrap_into_box(coords[atom][0], box_size[0])
    coords[atom][1] = wrap_into_box(coords[atom][1], box_size[1])
    coords[atom][2] = wrap_into_box(coords[atom][2], box_size[2])

    # 计算新能量
    new_energy = calculate_energy()

    accept = False

    # 如果能量降低，则自动接受本次移动
    if (new_energy <= old_energy):
        accept = True

    else:
        # 执行蒙特卡洛测试——比较
        # exp( -(E_new - E_old) / kT ) >= rand(0,1)
        x = math.exp( -(new_energy - old_energy) / kT )

        if (x >= random.uniform(0.0,1.0)):
            accept = True
        else:
            accept = False

    if accept:
        # 接受移动
        naccept += 1
        total_energy = new_energy
    else:
        # 拒绝移动——储存旧坐标
        nreject += 1

        # 储存旧构象
        coords = copy.deepcopy(old_coords)

        total_energy = old_energy

    # 每10次移动打印一次能量
    if move % 10 == 0:
        print("%s %s  %s  %s" % (move, total_energy, naccept, nreject))


    # 每100次移动打印一次坐标
    if move % 100 == 0:
        print_pdb(move)