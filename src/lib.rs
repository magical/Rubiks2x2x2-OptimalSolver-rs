#![allow(non_upper_case_globals)]
#![allow(non_snake_case)]
#![allow(dead_code)]
#![allow(unused_imports)]

mod coord {
    pub struct CoordCube {

    }

    impl CoordCube {
        pub fn new() -> CoordCube {
            return CoordCube {}
        }
    }
}

mod enums {
    /**

      The names of the facelet positions of the cube

                |********|
                |*U1**U2*|
                |********|
                |*U3**U4*|
                |********|
        |********|********|********|********|
        |*L1**L2*|*F1**F2*|*R1**R2*|*B1**B2*|
        |********|********|********|********|
        |*L3**L4*|*F3**F4*|*R3**R4*|*B3**B4*|
        |********|********|********|********|
                |********|
                |*D1**D2*|
                |********|
                |*D3**D4*|
                |********|

        A cube definition string "UBL..." means for example: In position U1 we have the U-color, in position U2 we have the
        B-color, in position U3 we have the L color etc. according to the order U1, U2, U3, U4, R1, R2, R3, R4, F1, F2, F3,
        F4, D1, D2, D3, D4, L1, L2, L3, L4, B1, B2, B3, B4 of the enum constants.
    */
    #[derive(Eq,PartialEq,Debug,Copy,Clone)]
    pub enum Facelet {
        U1, U2, U3, U4,
        R1, R2,R3, R4,
        F1, F2, F3, F4,
        D1, D2, D3, D4,
        L1, L2, L3, L4,
        B1, B2, B3, B4,
    }

    /** The possible colors of the cube facelets. Color U refers to the color of the U(p)-face etc.
        Also used to name the faces itself. **/
    #[derive(Eq,PartialEq,Debug,Copy,Clone)]
    pub enum Color {
        U, R, F, D, L, B,
    }

    /** The names of the corner positions of the cube. Corner URF e.g. has an U(p), a R(ight) and a F(ront) facelet. */
    #[derive(Eq,PartialEq,Debug,Copy,Clone)]
    pub enum Corner {
        URF, UFL, ULB, UBR, DRB, DFR, DLF, DBL,
    }

    use std::slice::Iter;
    impl Corner {
        pub fn iter() -> Iter<'static, Corner> {
            use self::Corner::*;
            static CORNERS: [Corner;8] = [URF, UFL, ULB, UBR, DRB, DFR, DLF, DBL];
            return CORNERS.into_iter()
        }
    }

    /** The moves in the faceturn metric. Not to be confused with the names of the facelet positions in class Facelet. */
    #[derive(Eq,PartialEq,Debug,Copy,Clone)]
    pub enum Move {
        U1, U2, U3, R1, R2, R3, F1, F2, F3
    }

}

mod defs {
    use enums::Facelet as Fc;
    use enums::Color as Cl;

    // Map the corner positions to facelet positions.
    pub static cornerFacelet: [[Fc;3];8] = [
    [Fc::U4, Fc::R1, Fc::F2], [Fc::U3, Fc::F1, Fc::L2], [Fc::U1, Fc::L1, Fc::B2], [Fc::U2, Fc::B1, Fc::R2],
                     [Fc::D4, Fc::R4, Fc::B3], [Fc::D2, Fc::F4, Fc::R3], [Fc::D1, Fc::L4, Fc::F3], [Fc::D3, Fc::B4, Fc::L3],
                     ];

    // Map the corner positions to facelet colors.
    pub static cornerColor: [[Cl;3];8]=[[Cl::U, Cl::R, Cl::F], [Cl::U, Cl::F, Cl::L], [Cl::U, Cl::L, Cl::B], [Cl::U, Cl::B, Cl::R],
                   [Cl::D, Cl::R, Cl::B], [Cl::D, Cl::F, Cl::R], [Cl::D, Cl::L, Cl::F], [Cl::D, Cl::B, Cl::L]
                   ];

    pub const N_MOVE: u32 = 9; // number of possible face moves
    pub const N_TWIST: u32 = 729; // 3^6 possible corner orientations
    pub const N_CORNERS: u32 = 5040; // 7! corner permutations in phase 2

}

mod cubie {
    //! The 2x2x2 cube on the cubie level is described by the permutation and orientations of the corners

    use defs::{cornerFacelet, cornerColor, N_CORNERS, N_TWIST};
    use enums::{Color,Corner as Co};
    use misc::{rotate_left,rotate_right};
    extern crate rand;

    fn randrange(a: u32, b: u32) -> u32 {
        use cubie::rand::Rng;
        let num = rand::thread_rng().gen_range(a, b);
        return num;
    }

    // the basic six cube moves described by permutations and changes in orientation

    // Up-move
    static cpU: [Co;8] = [Co::UBR, Co::URF, Co::UFL, Co::ULB,  Co::DRB, Co::DFR, Co::DLF, Co::DBL];
    static coU: [i32;8] = [0, 0, 0, 0, 0, 0, 0, 0];

    // Right-move
    static cpR: [Co;8] = [Co::DFR, Co::UFL, Co::ULB, Co::URF, Co::UBR, Co::DRB, Co::DLF, Co::DBL];  // permutation of the corners
    static coR: [i32;8] = [2, 0, 0, 1, 2, 1, 0, 0];  // changes of the orientations of the corners

    // Front-move
    static cpF: [Co;8] = [Co::UFL, Co::DLF, Co::ULB, Co::UBR, Co::DRB, Co::URF, Co::DFR, Co::DBL];
    static coF: [i32;8] = [1, 2, 0, 0, 0, 2, 1, 0];

    const CUBE_OK: bool = true;


    /// Represents a 2x2x2 cube on the cubie level with 8 corner cubies and the corner orientations.
    ///
    /// Is also used to represent the 18 cube moves.
    #[derive(Eq,PartialEq,Debug,Clone)]
    pub struct CubieCube {
        cp: [Co;8], // corner permutation
        co: [i32;8], // corner orientation
    }

    impl CubieCube {
        /// Initializes corners and edges.
        /// :param cp: corner permutation
        /// :param co: corner orientation
        pub fn new(cp: Option<[Co;8]>, co: Option<[i32;8]>) -> CubieCube {
            CubieCube{
                cp: if cp.is_none() {
                    [Co::URF, Co::UFL, Co::ULB, Co::UBR, Co::DRB, Co::DFR, Co::DLF, Co::DBL]
                } else {
                    cp.unwrap()
                },
                co: if co.is_none() {
                    [0; 8]
                } else {
                    co.unwrap()
                },
            }
        }

        /// Prints string for a cubie cube.
        pub fn to_str(&self) -> String {
            let mut s = String::new();
            for i in 0..self.cp.len() {
                s += &format!("({:?}, {})", self.cp[i], self.co[i]);
            }
            return s;
        }

        /*
        /// Returns a facelet representation of the cube.
        pub fn to_facelet_cube(&self) -> FaceCube {
            let fc = face::FaceCube()
            for i in 0..8 {
                let j = self.cp[i]  # corner j is at corner position i
                let ori = self.co[i]  # orientation of C j at position i
                for k in 0..3 {
                    fc.f[cornerFacelet[i][(k + ori) % 3]] = cornerColor[j][k]
                }
            }
            return fc
        }
        */

        /// Multiplies this cubie cube with another cubie cube b. Does not change b.
        pub fn multiply(&self, b: &Self) -> Self {
            let mut c_perm = [Co::URF; 8];
            let mut c_ori = [0; 8];
            let mut ori: i32 = 0;
            for c in 0..8 {
                c_perm[c] = self.cp[b.cp[c] as usize];
                let ori_a = self.co[b.cp[c] as usize];
                let ori_b = b.co[c];
                if ori_a < 3 && ori_b < 3 {  // two regular cubes
                    ori = ori_a + ori_b;
                    if ori >= 3 {
                        ori -= 3;
                    }
                }
                else if ori_a < 3 && 3 <= ori_b {  // cube b is in a mirrored state
                    ori = ori_a + ori_b;
                    if ori >= 6 {
                        ori -= 3;  // the composition also is in a mirrored state
                    }
                }
                else if ori_a >= 3 && 3 > ori_b {  // cube a is in a mirrored state
                    ori = ori_a - ori_b;
                    if ori < 3 {
                        ori += 3;  // the composition is a mirrored cube
                    }
                }
                else if ori_a >= 3 && ori_b >= 3 {  // if both cubes are in mirrored states
                    ori = ori_a - ori_b;
                    if ori < 0 {
                        ori += 3; // the composition is a regular cube
                    }
                }
                c_ori[c] = ori;
            }

            return Self{
                cp: c_perm,
                co: c_ori,
            }
        }

        /// Stores the inverse of this cubie cube in d.
        pub fn inv_cubie_cube(&self, d: &mut Self) {
            for cref in Co::iter() {
                let c: Co = *cref;
                d.cp[self.cp[c as usize] as usize] = c;
            }
            for c in 0..8 {
                let ori = self.co[d.cp[c] as usize];
                if ori >= 3 {
                    d.co[c] = ori;
                }
                else {
                    d.co[c] = -ori;
                    if d.co[c] < 0 {
                        d.co[c] += 3;
                    }
                }
            }
        }


        //
        // coordinates for 2x2x2 cube
        //

        /// The twist of the 8 corners. 0 <= twist < 729. The DBL-corner is fixed.
        pub fn get_corntwist(&self) -> u32 {
            let mut ret = 0;
            for i in 0..8 {
                ret = 3 * ret + (self.co[i] as u32);
            }
            return ret;
        }

        pub fn set_cornertwist(&mut self, mut twist: u32) {
            let mut twistparity = 0;
            for i in (0..8).rev() {
                self.co[i] = (twist % 3) as i32;
                twistparity += self.co[i];
                twist /= 3;
            }
            self.co[Co::DLF as usize] = (3 - twistparity % 3) % 3;
            // XXX need mathematical mod?
        }



        /// The permutation of the 8 corners. 0 <= corners < 5040. The DLB_corner is fixed.
        pub fn get_cornperm(&self) -> u32{
            let mut perm = self.cp;  // duplicate cp
            let mut b = 0;
            for jref in Co::iter().rev() {
                let j = *jref;
                let mut k = 0;
                while perm[j as usize] != j {
                    rotate_left(&mut perm, 0, j as usize);
                    k += 1;
                }
                b = (j as u32 + 1) * b + k;
            }
            return b;
        }

        pub fn set_corners(&mut self, mut idx: u32) {
            for j in Co::iter() {
                self.cp[*j as usize] = *j;
            }
            for jref in Co::iter() {
                let j = *jref as u32;
                let mut k = idx % (j+ 1);
                idx /= j;
                while k > 0 {
                    rotate_right(&mut self.cp, 0, j as usize);
                    k -= 1;
                }
            }
        }

        //
        // end coordinates for 2x2x2 cube
        //

        // other usefull functions

        /*
        def randomize(self):
            """Generates a random cube. The probability is the same for all possible states."""
            self.set_corners(randrange(N_CORNERS))
            self.set_cornertwist(randrange(N_TWIST))

        */

        /// Checks if cubiecube is valid
        pub fn verify(&self) -> Result<(), &'static str> {
            let mut corner_count = [0; 8];
            for c in self.cp.iter() {
                corner_count[*c as usize] += 1;
            }
            for count in corner_count.iter() {
                if *count != 1 {
                    return Err("Error: Some corners are undefined.");
                }
            }

            let mut s = 0;
            for i in self.co.iter() {
                s += *i;
            }
            if s % 3 != 0 {
                return Err("Error: Total corner twist is wrong.");
            }

            return Ok(());
        }

    }



    //
    // these cubes represent the basic cube moves
    //

    pub static basicMoveCube: [&CubieCube; 3] = [
        //&CubieCube{ cp: cpU, co: coU },
        //&CubieCube{ cp: cpR, co: coR },
        //&CubieCube{ cp: cpF, co: coF },

        &CubieCube{ cp: [Co::UBR, Co::URF, Co::UFL, Co::ULB,  Co::DRB, Co::DFR, Co::DLF, Co::DBL], co:  [0, 0, 0, 0, 0, 0, 0, 0]},
        &CubieCube{ cp: [Co::DFR, Co::UFL, Co::ULB, Co::URF, Co::UBR, Co::DRB, Co::DLF, Co::DBL], co: [2, 0, 0, 1, 2, 1, 0, 0]},
        &CubieCube{ cp: [Co::UFL, Co::DLF, Co::ULB, Co::UBR, Co::DRB, Co::URF, Co::DFR, Co::DBL], co: [1, 2, 0, 0, 0, 2, 1, 0] },
    ];


/*
    // these cubes represent the all 9 cube moves
    fn moveCube() -> [CubieCube; 9] {

    }

moveCube = [0] * 9
for c1 in [Color.U, Color.R, Color.F]:
    cc = CubieCube()
    for k1 in range(3):
        cc.multiply(basicMoveCube[c1])
        moveCube[3 * c1 + k1] = CubieCube(cc.cp, cc.co)

*/

}

mod misc {
    pub fn rotate_left<T>(_arr: &mut[T], _l: usize, _r: usize) {}
    pub fn rotate_right<T>(_arr: &mut [T], _l: usize, _r: usize) {}
}

mod moves {
    use cubie as cb;
    use enums;
    use defs::{N_TWIST, N_CORNERS, N_MOVE};

    pub fn init() {
        let mut a = cb::CubieCube::new(None, None);

        let mut corntwist_move = [0 as u32; (N_TWIST * N_MOVE) as usize];
        for i in 0..N_TWIST {
            a.set_cornertwist(i);
            for jref in [enums::Color::U, enums::Color::R, enums::Color::F].iter() {
                let j = *jref;
                for k in 0..3 {
                    a = a.multiply(cb::basicMoveCube[j as usize]);
                    corntwist_move[(N_MOVE * i + 3 * (j as u32) + k) as usize] = a.get_corntwist();
                }
                a = a.multiply(cb::basicMoveCube[j as usize]);
            }
        }

        // TODO: cache as file
        let mut cornperm_move = [0 as u32; (N_CORNERS * N_MOVE) as usize];
        for i in 0..N_CORNERS {
            if (i+1)%200 == 0 {
                print!(".")
            }
            a.set_corners(i);
            for jref in [enums::Color::U, enums::Color::R, enums::Color::F].iter() {
                let j = *jref;
                for k in 0..3 {
                    a = a.multiply(cb::basicMoveCube[j as usize]);
                    cornperm_move[(N_MOVE * i + 3 * (j as u32) + k) as usize] = a.get_cornperm();
                }
                a = a.multiply(cb::basicMoveCube[j as usize]);
            }
        }
    }
}
