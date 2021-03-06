#![allow(non_upper_case_globals)]
#![allow(non_snake_case)]
#![allow(dead_code)]

    use std::sync::{Once, ONCE_INIT};
    use std::slice::Iter;

    // some "constants"
    const N_MOVE: u32 = 9; // number of possible face moves
    const N_TWIST: u32 = 729; // 3^6 possible corner orientations
    const N_CORNERS: u32 = 5040; // 7! corner permutations in phase 2

    //
    // The cube on the coordinate level is described by a 3-tuple of natural numbers in phase 1 and phase 2.
    //

    const SOLVED: u32 = 0; // 0 is index of solved state (except for u_edges coordinate)

    /// Represents a 2x2x2 cube on the coordinate level.
    ///
    /// A state is uniquely determined by the coordinates corntwist and cornperm.
    #[derive(Eq,PartialEq,Debug,Clone)]
    struct CoordCube {
        corntwist: u32, // twist of corners
        cornperm: u32, // permutations of corners
    }

    impl CoordCube {
        fn new() -> CoordCube {
            return CoordCube { corntwist: SOLVED, cornperm: SOLVED }
        }

        fn from_cubie_cube(cc: &CubieCube) -> CoordCube {
            return CoordCube { corntwist: cc.get_corntwist(), cornperm: cc.get_cornperm() }
        }

        fn to_str(&self) -> String {
            return format!("(corntwist: {}, cornperm: {})", self.corntwist, self.cornperm)
        }
    }

    //
    // Enumerations which improve the readability of the code
    //

    /**

      The names of the facelet positions of the cube

      ```text

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

        ```

      A cube definition string "UBL..." means for example: In position U1 we have the U-color, in position U2 we have the
      B-color, in position U3 we have the L color etc. according to the order U1, U2, U3, U4, R1, R2, R3, R4, F1, F2, F3,
      F4, D1, D2, D3, D4, L1, L2, L3, L4, B1, B2, B3, B4 of the enum constants.
    */
    #[derive(Eq,PartialEq,Debug,Copy,Clone)]
    enum Facelet {
        U1, U2, U3, U4,
        R1, R2, R3, R4,
        F1, F2, F3, F4,
        D1, D2, D3, D4,
        L1, L2, L3, L4,
        B1, B2, B3, B4,
    }

    /** The possible colors of the cube facelets. Color U refers to the color of the U(p)-face etc.
        Also used to name the faces itself. **/
    #[derive(Eq,PartialEq,Debug,Copy,Clone)]
    enum Color {
        U, R, F, D, L, B,
    }

    impl Color {
        fn iter() -> Iter<'static, Color> {
            use Color::*;
            static COLORS: [Color;6] = [U, R, F, D, L, B];
            return COLORS.into_iter()

        }
    }

    /** The names of the corner positions of the cube. Corner URF e.g. has an U(p), a R(ight) and a F(ront) facelet. */
    #[derive(Eq,PartialEq,Debug,Copy,Clone)]
    enum Corner {
        URF, UFL, ULB, UBR, DRB, DFR, DLF, DBL,
    }

    impl Corner {
        fn iter() -> Iter<'static, Corner> {
            use Corner::*;
            static CORNERS: [Corner;8] = [URF, UFL, ULB, UBR, DRB, DFR, DLF, DBL];
            return CORNERS.into_iter()
        }
    }

    /** The moves in the faceturn metric. Not to be confused with the names of the facelet positions in class Facelet. */
    #[derive(Eq,PartialEq,Debug,Copy,Clone)]
    enum Move {
        U1, U2, U3, R1, R2, R3, F1, F2, F3
    }

    impl Move {
        fn iter() -> Iter<'static, Move> {
            use Move::*;
            static MOVES: [Move;9] = [U1, U2, U3, R1, R2, R3, F1, F2, F3];
            return MOVES.into_iter()
        }

        fn name(self) -> &'static str {
            use Move::*;
            return match self {
                U1=>"U", U2=> "U2", U3=>"U'",
                F1=>"F", F2=> "F2", F3=>"F'",
                R1=>"R", R2=> "R2", R3=>"R'",
            }
        }
    }

    //
    // some definitions
    //

    use self::Facelet as Fc;
    use self::Color as Cl;

    // Map the corner positions to facelet positions.
    static cornerFacelet: [[Facelet;3];8] = [
        [Fc::U4, Fc::R1, Fc::F2], [Fc::U3, Fc::F1, Fc::L2], [Fc::U1, Fc::L1, Fc::B2], [Fc::U2, Fc::B1, Fc::R2],
        [Fc::D4, Fc::R4, Fc::B3], [Fc::D2, Fc::F4, Fc::R3], [Fc::D1, Fc::L4, Fc::F3], [Fc::D3, Fc::B4, Fc::L3],
    ];

    // Map the corner positions to facelet colors.
    static cornerColor: [[Cl;3];8]=[
        [Cl::U, Cl::R, Cl::F], [Cl::U, Cl::F, Cl::L], [Cl::U, Cl::L, Cl::B], [Cl::U, Cl::B, Cl::R],
        [Cl::D, Cl::R, Cl::B], [Cl::D, Cl::F, Cl::R], [Cl::D, Cl::L, Cl::F], [Cl::D, Cl::B, Cl::L]
    ];

    //
    // The 2x2x2 cube on the facelet level is described by positions of the colored stickers.
    //

    static default_faces: [Color; 4*6] = [
        Color::U, Color::U, Color::U, Color::U,
        Color::R, Color::R, Color::R, Color::R,
        Color::F, Color::F, Color::F, Color::F,
        Color::D, Color::D, Color::D, Color::D,
        Color::L, Color::L, Color::L, Color::L,
        Color::B, Color::B, Color::B, Color::B,
    ];

    /// Represents a 2x2x2 cube on the facelet level with 24 colored facelets
    struct FaceCube {
        f: [Color; 4*6]
    }


    impl FaceCube {
        fn new(f: [Color; 4*6]) -> FaceCube {
            FaceCube{ f: f }
        }
        /// Constructs a facelet cube from a string. See class Facelet(IntEnum) in enums.py for string format.
        ///
        /// The color scheme is detected automatically.
        fn from_string(s: &str) -> Result<FaceCube, &'static str> {
            if s.len() < 24 {
                return Err("Error: Cube definition string contains less than 24 facelets")
            }
            if s.len() > 24 {
                return Err("Error: Cube definition string contains more than 24 facelets")
            }

            let mut f = default_faces.clone();
            let mut cnt = [0; 6];
            for (i, ch) in s.chars().enumerate() {
                if ch == 'U' {
                    f[i] = Color::U;
                    cnt[Color::U as usize] += 1;
                } else if ch == 'R' {
                    f[i] = Color::R;
                    cnt[Color::R as usize] += 1;
                } else if ch == 'F' {
                    f[i] = Color::F;
                    cnt[Color::F as usize] += 1;
                } else if ch == 'D' {
                    f[i] = Color::D;
                    cnt[Color::D as usize] += 1;
                } else if ch == 'L' {
                    f[i] = Color::L;
                    cnt[Color::L as usize] += 1;
                } else if ch == 'B' {
                    f[i] = Color::B;
                    cnt[Color::B as usize] += 1;
                }
            }
            if !cnt.iter().all(|x| *x == 4) {
                return Err("Error: Cube definition string does not contain exactly 4 facelets of each color.")
            }

            // remap colors if necessary

            let mut col = [Color::U; 3];
            for i in 0..3 {
                col[i] = f[cornerFacelet[Corner::DBL as usize][i] as usize]  // colors of the DBL-corner
            }
            let mut map_col = [None; 6];
            for i in 0..3 {
                map_col[col[i] as usize] = Some(cornerColor[Corner::DBL as usize][i])  // map colors to right colors
            }
            // now remap the remaining colors, try all possibilites
            let a = vec![
                [Color::U, Color::R, Color::F], [Color::U, Color::F, Color::R], [Color::R, Color::U, Color::F],
                [Color::R, Color::F, Color::U], [Color::F, Color::U, Color::R], [Color::F, Color::R, Color::U]];
            let mut empty = vec![];
            for i in Color::iter() {
                let i = *i;
                if map_col[i as usize].is_none() {
                    empty.push(i as usize)  // empty contains the 3 indices of the yet nonmapped colors
                }
            }

            let fsave = f;
            for c in a {
                for i in 0..3 {
                    map_col[empty[i]] = Some(c[i]);
                }
                let mut f = fsave.clone();
                for i in 0..24 {
                    f[i] = map_col[fsave[i] as usize].unwrap();  // remap the colors
                }
                let fc = FaceCube{f: f};
                let cc = fc.to_cubie_cube();
                if cc.verify().is_ok() {
                    return Ok(fc);
                }
            }

            return Err("Error: Facelet configuration does not define a valid cube.")

        }

        /// Gives string representation of the facelet cube.
        fn to_string(&self) -> String {
            self.f.iter().map(|c|
                match *c {
                    Color::U => "U",
                    Color::R => "R",
                    Color::F => "F",
                    Color::D => "D",
                    Color::L => "L",
                    Color::B => "B",
                }
            ).collect()
        }

        /// Returns a cubie representation of the facelet cube.
        fn to_cubie_cube(&self) -> CubieCube {
            let mut cp = [Corner::URF; 8]; // invalidate corner permutation
            let mut co = [0 as i32; 8];

            for i in Corner::iter() {
                let i = *i;
                let mut col0 = None;
                let mut ori0: i32 = 2;
                let fac = cornerFacelet[i as usize]; // facelets of corner  at position i
                for ori in 0..3 {
                    if self.f[fac[ori as usize] as usize] == Color::U || self.f[fac[ori as usize] as usize] == Color::D {
                        col0 = Some(self.f[fac[ori as usize] as usize]);
                        ori0 = ori;
                        break
                    }
                }
                let col1 = self.f[fac[((ori0+1)%3) as usize] as usize]; // colors which identify the corner at position i
                let col2 = self.f[fac[((ori0+2)%3) as usize] as usize];
                for j in Corner::iter() {
                    let j = *j;
                    let col = cornerColor[j as usize]; // colors of corner j
                    if col0 == Some(col[0]) && col1 == col[1] && col2 == col[2]  {
                        cp[i as usize] = j; // we have coner j in corner position i
                        co[i as usize] = ori0;
                    }
                }
            }

            return CubieCube{ cp: cp, co: co }
        }
    }


    //
    // The 2x2x2 cube on the cubie level is described by the permutation and orientations of the corners
    //

    use self::Corner as Co;


    /// Represents a 2x2x2 cube on the cubie level with 8 corner cubies and the corner orientations.
    ///
    /// Is also used to represent the 18 cube moves.
    #[derive(Eq,PartialEq,Debug,Clone)]
    struct CubieCube {
        cp: [Corner;8], // corner permutation
        co: [i32;8], // corner orientation
    }

    impl CubieCube {
        /// Constructs a cube in the solved position.
        fn new() -> CubieCube {
            CubieCube{
                cp: [Co::URF, Co::UFL, Co::ULB, Co::UBR, Co::DRB, Co::DFR, Co::DLF, Co::DBL],
                co: [0; 8],
            }
        }

        /// Prints string for a cubie cube.
        fn to_str(&self) -> String {
            let mut s = String::new();
            for i in 0..self.cp.len() {
                s += &format!("({:?}, {})", self.cp[i], self.co[i]);
            }
            return s;
        }

        /// Returns a facelet representation of the cube.
        fn to_facelet_cube(&self) -> FaceCube {
            let mut f = [Color::U; 24];
            for i in 0..8 {
                let j = self.cp[i];  // corner j is at corner position i
                let ori = self.co[i];  // orientation of C j at position i
                for k in 0..3 {
                    f[cornerFacelet[i][((k + ori) % 3) as usize] as usize] = cornerColor[j as usize][k as usize]
                }
            }
            return FaceCube::new(f);
        }

        /// Multiplies this cubie cube with another cubie cube b. Does not change b.
        fn multiply(&self, b: &Self) -> Self {
            let mut c_perm = [Corner::URF; 8];
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

        /// Returns inverse of this cubie cube.
        fn inv(&self) -> Self {
            let mut d = Self::new();
            for c in Corner::iter() {
                let c = *c;
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
            return d;
        }


        //
        // coordinates for 2x2x2 cube
        //

        /// The twist of the 8 corners. 0 <= twist < 729. The DBL-corner is fixed.
        fn get_corntwist(&self) -> u32 {
            let mut ret = 0;
            for i in 0..6 { // note: last two corners not included
                ret = 3 * ret + (self.co[i] as u32);
            }
            assert!(ret < 729);
            return ret;
        }

        fn set_cornertwist(&mut self, mut twist: u32) {
            assert!(twist < 729);
            let mut twistparity = 0;
            for i in (0..6).rev() {
                self.co[i] = (twist % 3) as i32;
                twistparity += self.co[i];
                twist /= 3;
            }
            self.co[Corner::DLF as usize] = (3 - twistparity % 3) % 3;
            // XXX need mathematical mod?
        }

        /// The permutation of the 8 corners. 0 <= corners < 5040. The DLB_corner is fixed.
        fn get_cornperm(&self) -> u32{
            let mut perm = self.cp;  // duplicate cp
            let mut b = 0;
            for j in Corner::iter().rev() {
                let j = *j;
                let mut k = 0;
                while perm[j as usize] != j {
                    rotate_left(&mut perm[0..(j as usize + 1)]);
                    k += 1;
                }
                b = (j as u32 + 1) * b + k;
            }
            assert!(b < 5040);
            return b;
        }

        fn set_corners(&mut self, mut idx: u32) {
            assert!(idx < 5040);
            for j in Corner::iter() {
                self.cp[*j as usize] = *j;
            }
            for j in Corner::iter() {
                let j = *j as u32;
                let mut k = idx % (j+ 1);
                idx /= j + 1;
                while k > 0 {
                    rotate_right(&mut self.cp[0..(j as usize + 1)]);
                    k -= 1;
                }
            }
        }

        //
        // end coordinates for 2x2x2 cube
        //

        // other useful functions

        /*
        def randomize(self):
            """Generates a random cube. The probability is the same for all possible states."""
            self.set_corners(randrange(N_CORNERS))
            self.set_cornertwist(randrange(N_TWIST))

        */

        /// Checks if cubiecube is valid
        fn verify(&self) -> Result<(), &'static str> {
            let mut corner_count = [0; 8];
            for c in self.cp.iter() {
                corner_count[*c as usize] += 1;
            }
            for count in corner_count.iter() {
                if *count != 1 {
                    return Err("Error: Some corners are undefined.");
                }
            }

            let s: i32 = self.co.iter().sum();
            if s % 3 != 0 {
                return Err("Error: Total corner twist is wrong.");
            }

            return Ok(());
        }

    }

    /// Rotates array arr right.
    fn rotate_left<T>(arr: &mut[T]) where T: Copy {
        let r = arr.len() - 1;
        let temp = arr[r];
        for i in (0..r).rev() {
            arr[i+1] = arr[i]
        }
        arr[0] = temp
    }

    /// Rotates array arr left
    fn rotate_right<T>(arr: &mut [T]) where T: Copy {
        let r = arr.len() - 1;
        let temp = arr[0];
        for i in 0..r {
            arr[i] = arr[i+1]
        }
        arr[r] = temp
    }

    //
    // these cubes represent the basic cube moves
    //

    // the basic six cube moves described by permutations and changes in orientation

    static basicMoveCube: [&CubieCube; 3] = [
        // Up-move
        &CubieCube{ cp: [Co::UBR, Co::URF, Co::UFL, Co::ULB, Co::DRB, Co::DFR, Co::DLF, Co::DBL], co:  [0, 0, 0, 0, 0, 0, 0, 0] },

        // Right-move
        &CubieCube{ cp: [Co::DFR, Co::UFL, Co::ULB, Co::URF, Co::UBR, Co::DRB, Co::DLF, Co::DBL], co: [2, 0, 0, 1, 2, 1, 0, 0] },

        // Front-move
        &CubieCube{ cp: [Co::UFL, Co::DLF, Co::ULB, Co::UBR, Co::DRB, Co::URF, Co::DFR, Co::DBL], co: [1, 2, 0, 0, 0, 2, 1, 0] },
    ];


    //
    // Movetables describe the transformation of the coordinates by cube moves.
    //

    // Move table for the the corners.

    /// Group performs operations in the Rubiks cube group.
    struct Group {
        corntwist_move: [u16; n_corntwist],
        cornperm_move: [u16; n_cornperm],
    }

    impl Group {
        fn apply_move(&self, c: &CoordCube, m: Move) -> CoordCube {
            CoordCube{
                cornperm: self.cornperm_move[(N_MOVE*c.cornperm + m as u32) as usize] as u32,
                corntwist: self.corntwist_move[(N_MOVE*c.corntwist + m as u32) as usize] as u32,
            }
        }
    }

    const n_corntwist: usize = (N_TWIST * N_MOVE) as usize;
    const n_cornperm: usize = (N_CORNERS*N_MOVE) as usize;

    static mut GROUP: Group = Group{
        // The twist coordinate describes the 3^6 = 729 possible orientations of the 8 corners
        corntwist_move: [0; n_corntwist],

        // The corners coordinate describes the 7! = 5040 permutations of the corners.
        cornperm_move: [0; n_cornperm],
    };

    fn get_group() -> &'static Group {
        static once: Once = ONCE_INIT;
        once.call_once(|| {
            init_corntwist(unsafe { &mut GROUP.corntwist_move });
            init_cornperm(unsafe { &mut GROUP.cornperm_move });
        });
        return unsafe{ &GROUP };
    }

    fn init_corntwist(corntwist_move: &mut [u16; n_corntwist]) {
        println!("creating move_corntwist table");
        let mut a = CubieCube::new();
        for i in 0..N_TWIST {
            a.set_cornertwist(i);
            for j in [Color::U, Color::R, Color::F].iter() { // three faces U, R, F
                let j = *j;
                for k in 0..3 { // three moves for each face, for example U, U2, U3 = U'
                    a = a.multiply(basicMoveCube[j as usize]);
                    let idx = (N_MOVE*(i as u32) + 3*(j as u32) + k) as usize;
                    corntwist_move[idx] = a.get_corntwist() as u16;
                }
                a = a.multiply(basicMoveCube[j as usize]); // 4. move restores face
            }
        }
    }

    fn init_cornperm(cornperm_move: &mut [u16; (N_CORNERS*N_MOVE) as usize]) {
        println!("creating move_cornperm table");
        let mut a = CubieCube::new();
        // TODO: cache as file
        // Move table for the corners. corner < 40320
        for i in 0..N_CORNERS {
            a.set_corners(i);
            for j in [Color::U, Color::R, Color::F].iter() { // three faces U, R, F
                let j = *j;
                for k in 0..3 {
                    a = a.multiply(basicMoveCube[j as usize]);
                    let idx = (N_MOVE*(i as u32) + 3*(j as u32) + k) as usize;
                    cornperm_move[idx] = a.get_cornperm() as u16;
                }
                a = a.multiply(basicMoveCube[j as usize]);
            }
        }
    }

    //
    // The pruning table cuts the search tree during the search.
    // In this case it it gives the exact distance to the solved state.
    //

    const n_cornerprun: usize = (N_CORNERS * N_TWIST) as usize;
    static mut cornerprun_table: [i8; n_cornerprun] = [-1; n_cornerprun];

    fn get_pruning_table() -> &'static [i8; n_cornerprun] {
        static once: Once = ONCE_INIT;
        unsafe {
            once.call_once(|| { init_cornerprun_table(&mut cornerprun_table); });
            return &cornerprun_table;
        }
    }

    /// creates/loads the corner pruning table
    fn init_cornerprun_table(corner_depth: &mut [i8; n_cornerprun]) {

        let group = get_group();

        println!("creating cornerprun table...");
        corner_depth[0] = 0;

        let mut done = 1;
        let mut depth = 0;
        while done != N_CORNERS * N_TWIST {
            for corners in 0..N_CORNERS {
                for twist in 0..N_TWIST {
                    let cube = CoordCube{ cornperm: corners, corntwist: twist };
                    if corner_depth[(N_TWIST * corners + twist) as usize] == depth {
                        for m in Move::iter() { // Move
                            let cube1 = group.apply_move(&cube, *m);
                            let idx1 = (N_TWIST * cube1.cornperm + cube1.corntwist) as usize;
                            if corner_depth[idx1] == -1 { // entry not yet filled
                                corner_depth[idx1] = depth+1;
                                done += 1;
                            }
                        }
                    }
                }
            }
            depth += 1
        }
    }

    //
    // The solve function computes all optimal solving maneuvers
    //

    struct Solver<'a> {
        solutions: Vec<Vec<Move>>,
        sofar: Vec<Move>,
        group: &'a Group,
        corner_depth: &'a [i8; n_cornerprun],
    }

    impl<'a> Solver<'a> {
        fn search(&mut self, cube: CoordCube, togo: i8) {
            if togo == 0 {
                self.solutions.push(self.sofar.clone())
            } else {
                for m in Move::iter() {
                    let m = *m;
                    if self.sofar.len() > 0 {
                        if self.sofar[self.sofar.len()-1] as u32 / 3 == (m as u32) / 3 { // successive moves on same face
                            continue
                        }
                    }

                    let cube_new = self.group.apply_move(&cube, m);
                    if self.corner_depth[(N_TWIST * cube_new.cornperm + cube_new.corntwist) as usize] >= togo {
                        continue // impossible to reach solved cube in togo - 1 moves
                    }

                    self.sofar.push(m);
                    self.search(cube_new, togo - 1);
                    self.sofar.pop();
                }
            }
        }
    }

    /// Solves a 2x2x2 cube defined by its cube definition string.
    ///
    /// :param cubestring: The format of the string is given in the Facelet class defined in the file enums.py
    /// :return A list of all optimal solving maneuvers
    pub fn solve(cubestring: &str) -> Result<String,&'static str> {
        let fc = try!(FaceCube::from_string(cubestring));
        let cc = fc.to_cubie_cube();
        let _ = try!(cc.verify());

        let mut solver = Solver {
            solutions: Vec::new(),
            sofar: Vec::new(),

            group: get_group(),
            corner_depth: get_pruning_table(),
        };

        let co_cube = CoordCube::from_cubie_cube(&cc);
        let togo = solver.corner_depth[(N_TWIST * co_cube.cornperm + co_cube.corntwist) as usize];
        solver.search(co_cube, togo);

        let solutions = solver.solutions;

        let mut s = String::new();
        for sol in solutions.iter() {
            let mut ps = String::new();
            for m in sol {
                ps += m.name();
                ps += " ";
            }
            ps += &format!("({}f)\n", sol.len());
            s += &ps;
        }
        return Ok(s);
    }

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        println!("starting");
        let _ = super::solve("UUUURRRRFFFFDDDDLLLLBBBB");
        println!("ok");
    }
}
