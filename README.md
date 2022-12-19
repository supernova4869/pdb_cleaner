# pdb_cleaner
Made by Supernova (https://www.zhihu.com/people/zhang-jia-xing-42-34)

## 功能:
由于 DS 的默认蛋白原子命名规则与 gromacs 的不同, 通过 Prepare Protein 功能得到的 pdb 文件, gromacs 是不能直接处理的.
这里写了一个小脚本, 将 DS prepare 得到的 protein 修改为 gromacs 可接受的形式. 具体修改内容如下: 

1. 重命名His, 根据H的类型将His重命名为HID, HIE and HIP
2. 重命名Asp, 若为质子化状态(包含HD2)则重命名为ASH
3. 重命名Glu, 若为质子化状态(包含HE2)则重命名为GLH
4. 重命名末端羧基O, 将DS输出的OCT[12]改为OC[12]

## Usage
- Usage 0: pdb_cleaner # 直接查看此帮助
- Usage 1: pdb_cleaner input.pdb output.pdb
- Usage 2: pdb_cleaner input.pdb # 默认输出到new.pdb

## Examples

```bash
# get test.pdb file prepared by DS
pdb_cleaner test.pdb # generates new.pdb
gmx pdb2gmx -f new.pdb
```
