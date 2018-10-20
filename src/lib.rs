
#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}

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
    pub enum Facelet {
        U1, U2, U3, U4,
        R1, R2,R3, R4,
        F1, F2, F3, F4,
        D1, D2, D3, D4,
        L1, L2, L3, L4,
        B1, B2, B3, B4,
    }

    /** The names of the corner positions of the cube. Corner URF e.g. has an U(p), a R(ight) and a F(ront) facelet. */
    pub enum Color {
        U, R, F, D, L, B,
    }

    /** The moves in the faceturn metric. Not to be confused with the names of the facelet positions in class Facelet. */
    pub enum Move {}

}

mod defs {
    use enums::Facelet as Fc;
    use enums::Color as Cl;

    static cornerFacelet: [[Fc;3];8] = [
    [Fc::U4, Fc::R1, Fc::F2], [Fc::U3, Fc::F1, Fc::L2], [Fc::U1, Fc::L1, Fc::B2], [Fc::U2, Fc::B1, Fc::R2],
                     [Fc::D4, Fc::R4, Fc::B3], [Fc::D2, Fc::F4, Fc::R3], [Fc::D1, Fc::L4, Fc::F3], [Fc::D3, Fc::B4, Fc::L3],
                     ];

    static cornerColor: [[Cl;3];8]=[[Cl::U, Cl::R, Cl::F], [Cl::U, Cl::F, Cl::L], [Cl::U, Cl::L, Cl::B], [Cl::U, Cl::B, Cl::R],
                   [Cl::D, Cl::R, Cl::B], [Cl::D, Cl::F, Cl::R], [Cl::D, Cl::L, Cl::F], [Cl::D, Cl::B, Cl::L]
                   ];

    const N_MOVE: u32 = 9; // number of possible face moves
    const N_TWIST: u32 = 729; // 3^6 possible corner orientations
    const N_CORNERS: u32 = 5040; // 7! corner permutations in phase 2
}
