use std::{env, fmt, fs};
use std::fmt::Display;
use std::fs::File;
use std::io::{stdin, Write};

struct PDB {
    title: String,
    crystallographic: String,
    residues: Vec<Residue>,
}

struct Residue {
    atoms: Vec<Atom>,
}

impl Display for PDB {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut s = format!("REMARK {}\nCRYST1 {}\n", self.title, self.crystallographic);
        for res in &self.residues {
            for atom in &res.atoms {
                s.push_str(format!("{}", atom).as_str());
            }
        }
        write!(f, "{}", s)
    }
}

#[derive(Clone)]
enum AtomType {
    ATOM,
    HETATM,
    TER,
}


impl Display for AtomType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            AtomType::ATOM => write!(f, "ATOM  "),
            AtomType::HETATM => write!(f, "HETATM"),
            AtomType::TER => write!(f, "TER   ")
        }
    }
}

// contains TER
#[derive(Clone)]
struct Atom {
    typ: AtomType,
    atnum: i32,
    atid: String,
    resname: String,
    chainid: String,
    resid: i32,
    x: f64,
    y: f64,
    z: f64,
    occupy: f64,
    bf: f64,
    atname: String,
    charge: String,
}

impl Atom {
    fn new(typ: AtomType,
           atnum: i32,
           atid: String,
           resname: String,
           chainid: String,
           resid: i32,
           x: f64,
           y: f64,
           z: f64,
           occupy: f64,
           bf: f64,
           atname: String,
           charge: String) -> Atom {
        let mut atname = atname;
        if atname.is_empty() {
            atname = atid[0..1].to_string();
        }
        Atom { typ, atnum, atid, resname, chainid, resid, x, y, z, occupy, bf, atname, charge }
    }
}

impl Display for Atom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.typ {
            AtomType::TER => write!(f, "TER\n"),
            _ => write!(f, "{:<6}{:>5} {:4} {:3} {:1}{:4}    {:8.3}{:8.3}{:8.3}{:6.2}{:6.2} \
                         {:>12}{:<2}\n", self.typ, self.atnum, self.atid, self.resname, self.chainid,
                        self.resid, self.x, self.y, self.z, self.occupy, self.bf, self.atname, self.charge)
        }
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    match args.len() {
        1 => {
            println!("pdb_cleaner (by Supernova, https://www.zhihu.com/people/zhang-jia-xing-42-34)\n\n\
            功能:\n\
            将 DS prepare 得到的 protein 修改为 gromacs 可接受的形式\n\
            1. 重命名His, 根据H的类型将His重命名为HID, HIE and HIP\n\
            2. 重命名Asp, 若为质子化状态(包含HD2)则重命名为ASH\n\
            3. 重命名Glu, 若为质子化状态(包含HE2)则重命名为GLH\n\
            4. 重命名末端羧基O, 将DS输出的OCT[12]改为OC[12]\n\n\
            用法0: pdb_cleaner                        # 直接查看此帮助\n\
            用法1: pdb_cleaner input.pdb output.pdb\n\
            用法2: pdb_cleaner input.pdb              # 默认输出到new.pdb\n\n\
            按任意键退出...");
            stdin().read_line(&mut String::new()).unwrap();
        }
        2 => {
            let pdb = args[1].trim();
            clean_pdb(pdb, "new.pdb");
        }
        3 => {
            let pdb = args[1].trim();
            clean_pdb(pdb, &args[2]);
        }
        _ => {
            println!("错误输入, 直接执行pdb_cleaner看说明.\n\
            按任意键退出...");
            stdin().read_line(&mut String::new()).unwrap();
        }
    }
}

fn clean_pdb(pdb: &str, out: &str) {
    if let Ok(pdb) = fs::read_to_string(pdb) {
        let mut pdb = parse_pdb(pdb);
        fix_residues(&mut pdb);
        let mut out = File::create(out).unwrap();
        out.write_all(format!("{}", pdb).as_bytes()).unwrap();
    } else {
        println!("你这 pdb 有问题啊");
        stdin().read_line(&mut String::new()).unwrap();
    }
}

fn parse_pdb(pdb_content: String) -> PDB {
    let pdb_content: Vec<&str> = pdb_content.split("\n").collect();

    let mut title: String = pdb_content[0][7..].trim().to_string();
    if !title.ends_with("by supernova") {
        title.push_str(" by supernova");
    }
    let mut crystallographic = String::new();
    let mut residues: Vec<Residue> = vec![];
    let mut residue: Vec<Atom> = vec![];

    let mut cur_res: i32 = pdb_content[3][22..26].trim().parse().unwrap();
    for line in pdb_content {
        if line.starts_with("CRYST1") {
            crystallographic = line.trim().to_string();
        }
        if line.starts_with("ATOM") || line.starts_with("HETATM") || line.starts_with("TER") {
            let typ = match &line[..4] {
                "ATOM" => AtomType::ATOM,
                "HETA" => AtomType::HETATM,
                _ => AtomType::TER,
            };

            let mut atnum: i32 = 0;
            let mut atid: String = String::new();
            let mut resname: String = String::new();
            let mut chainid: String = String::new();
            let mut resid: i32 = 0;
            let mut x: f64 = 0.0;
            let mut y: f64 = 0.0;
            let mut z: f64 = 0.0;
            let mut occupy: f64 = 0.0;
            let mut bf: f64 = 0.0;
            let mut atname = String::new();
            let mut charge = String::new();
            if !line.starts_with("TER") {
                atnum = line[6..11].trim().parse().unwrap();
                atid = line[12..16].to_string();
                resname = line[17..20].trim().to_string();
                chainid = match line[21..22].trim().len() {
                    0 => "X".to_string(),
                    _ => line[21..22].trim().to_string()
                };
                resid = line[22..26].trim().parse().unwrap();
                x = line[30..38].trim().parse().unwrap();
                y = line[38..46].trim().parse().unwrap();
                z = line[46..54].trim().parse().unwrap();
                occupy = line[54..60].trim().parse().unwrap();
                bf = line[60..66].trim().parse().unwrap();
                if line.len() >= 77 {
                    atname.push_str(line[76..78].trim());
                    charge.push_str(line[78..80].trim());
                }
            }

            if line.starts_with("TER") {         // 到链终止
                residues.push(Residue { atoms: residue.to_vec() });
                residue.clear();
                residues.push(Residue {
                    atoms: vec![Atom {
                        typ: AtomType::TER,
                        atnum: 0,
                        atid: String::new(),
                        resname: String::new(),
                        chainid: String::new(),
                        resid: 0,
                        x: 0.0,
                        y: 0.0,
                        z: 0.0,
                        occupy: 0.0,
                        bf: 0.0,
                        atname: String::new(),
                        charge: String::new(),
                    }]
                });
            } else if resid == cur_res && !line.starts_with("END") {
                residue.push(Atom::new(
                    typ,
                    atnum,
                    atid,
                    resname,
                    chainid,
                    resid,
                    x,
                    y,
                    z,
                    occupy,
                    bf,
                    atname,
                    charge,
                ));
            } else {                                        // 到下一个残基
                cur_res = resid;
                residues.push(Residue { atoms: residue.to_vec() });
                residue.clear();
                residue.push(Atom::new(
                    typ,
                    atnum,
                    atid,
                    resname,
                    chainid,
                    resid,
                    x,
                    y,
                    z,
                    occupy,
                    bf,
                    atname,
                    charge,
                ));
            }
        };
    }
    PDB {
        title,
        crystallographic,
        residues,
    }
}

enum HisType {
    HID,
    HIE,
    HIP,
}

fn fix_residues(pdb: &mut PDB) {
    for res in &mut pdb.residues {
        if res.atoms.is_empty() {
            continue;
        }
        let is_his = match res.atoms[0].resname.as_str() {
            "HIS" => true,
            _ => false
        };
        let is_glu = match res.atoms[0].resname.as_str() {
            "GLU" => true,
            _ => false
        };
        let is_asp = match res.atoms[0].resname.as_str() {
            "ASP" => true,
            _ => false
        };

        if is_his {
            let mut ats: Vec<String> = vec![];
            for atom in &mut res.atoms {
                ats.push(atom.atid.to_string());
            }
            let his_type: HisType;
            if ats.contains(&" HE2".to_string()) && ats.contains(&" HD1".to_string()) {
                his_type = HisType::HIP;
            } else if ats.contains(&" HE2".to_string()) {
                his_type = HisType::HIE;
            } else {
                his_type = HisType::HID;
            }
            match his_type {
                HisType::HIP => for atom in &mut res.atoms {
                    atom.resname = "HIP".to_string();
                }
                HisType::HIE => for atom in &mut res.atoms {
                    atom.resname = "HIE".to_string();
                }
                HisType::HID => for atom in &mut res.atoms {
                    atom.resname = "HID".to_string();
                }
            }
        }

        if is_glu {
            let mut has_he2 = false;
            for atom in &mut res.atoms {
                if atom.atid == " HE2" {
                    has_he2 = true;
                }
            }
            if has_he2 {
                for atom in &mut res.atoms {
                    atom.resname = "GLH".to_string();
                }
            }
        }

        if is_asp {
            let mut has_hd2 = false;
            for atom in &res.atoms {
                if atom.atid == " HD2" {
                    has_hd2 = true;
                }
            }
            if has_hd2 {
                for atom in &mut res.atoms {
                    atom.resname = "ASH".to_string();
                }
            }
        }

        for atom in &mut res.atoms {
            // 末端羧基O命名
            if atom.atid == "1OCT" {
                atom.atid = " OC1".to_string();
            } else if atom.atid == "2OCT" {
                atom.atid = " OC2".to_string();
            }
        }
    }
}
