const V1: [usize; 4] = [0, 1, 2, 3];
const V2: [usize; 4] = [4, 5, 6, 7];
const ARROWS: usize = 8;
const VERTEX_NUM: usize = 2;
const VERTEX_MAP: [&[usize]; VERTEX_NUM] = [&V1, &V2];
const MAX_GENUS: usize = 6;


struct PermutationStructure {
	out: [usize; ARROWS],
}

fn is_in_array(item: usize, array: &[usize]) -> bool {

	for i in array {
		if *i==item {
			return true;
		}
	}

	return false;
}


fn check_connected(edge_map: &PermutationStructure) -> bool {

	// vertex keeps track of which vertices we have visited; a 0 means no visit
	// a 1 means visited but not checked, a 2 means checked.
	let mut vertex = [0; VERTEX_NUM];
	vertex[0] = 1;

	loop {
		// Find a vertex that has been visited and not checked.
		let mut v: usize = 0;
		while !(vertex[v]==1) {
			v += 1;
			if v == vertex.len() {
				break;
			}
		}

		if v < vertex.len() {
			// Mark vertex v as checked
			vertex[v] = 2;

			for idx in 0..VERTEX_MAP[v].len() {
				let a = (*VERTEX_MAP[v])[idx];
				let a2 = edge_map.out[a];

				let mut v2 = 0;
				loop {
					if is_in_array(a2, VERTEX_MAP[v2]) {
						break;
					} else {
						v2 += 1;
					}

				}

				if vertex[v2] == 0 {
					vertex[v2] = 1;
				}
			} 

		} else {
			break;
		}

	}


	if is_in_array(0, &vertex) {
		return false;
	} else {
		return true;
	}

}

fn convert_vertex_map() -> PermutationStructure {
	let mut ans = PermutationStructure	{
		out: [0; ARROWS],
	};

	for i in 0..ARROWS {
		let mut vertex = 0;
		loop { 
			if is_in_array(i, VERTEX_MAP[vertex]) {
				let mut a = 0;
				loop {
					if (*VERTEX_MAP[vertex])[a] == i {
						if a < VERTEX_MAP[vertex].len()-1 {
							ans.out[i] = (*VERTEX_MAP[vertex])[a+1];
						} else {
							ans.out[i] = (*VERTEX_MAP[vertex])[0];
						}
						break;
					} else {
						a += 1;
					}

				}
				break;

			} else {
				vertex += 1;
			}

		}
	}

	return ans;

}




fn product_permute(edge_map: &PermutationStructure, vertex_map: &PermutationStructure) -> PermutationStructure {
	let mut ans = PermutationStructure {
		out: [0; ARROWS],
	};

	for p in 0..ARROWS {
		ans.out[p] = vertex_map.out[edge_map.out[p]];
	}

	return ans;
}

fn face_count(permute: PermutationStructure) -> usize {
	let mut ans = 0;

	// travel represents each arrow in the graph;
	// we will follow each arrow until we complete a cycle;
	// then increment ans and move on to the next unused arrow.

	let mut travel = [false; ARROWS];

	loop {
		let mut a: usize = 0;
		while travel[a] {
			a += 1;
			if a == travel.len() {
				break;
			}
		}

		if a < travel.len() {
			let mut p = a;

			while !(travel[p]) {
				travel[p] = true;
				p = permute.out[p];

				if p == a {
					ans += 1;
					break;
				} 
			}
		} else {
			break;
		}




	}

	return ans;

}

fn count(genus: &mut [usize], edge_map: &PermutationStructure, vertex_map: &PermutationStructure) {
	if check_connected(edge_map) {
		let faces: isize = face_count(product_permute(edge_map, vertex_map)) as isize;
		let edges: isize = (ARROWS / 2) as isize;
		let vertices: isize = VERTEX_NUM as isize;

		let chi: isize = vertices - edges + faces;
		let idx: usize = ((2-chi)/2) as usize;
		genus[idx] += 1;
	}
}

// -------------------------------------------------------------------------
// Now for the real meat of the program. The following functions set up a nesting
// procedure for cycling through all possible choices for edge_map. 

fn nest2(nest_count: usize, which: &PermutationStructure, chooser: &mut PermutationStructure, 
	edge_map: &mut PermutationStructure, vertex_map: & PermutationStructure, 
	genus: &mut [usize] ) {

	if nest_count < ARROWS/2 - 1 {
		let idx0 = chooser.out[0];
		let which_idx = which.out[nest_count+1];
		let idx1 = chooser.out[which_idx];
		edge_map.out[idx0] = idx1;
		edge_map.out[idx1] = idx0;

		for j in 0..(ARROWS - 2*nest_count-1) {
			chooser.out[j] = chooser.out[j+1];
		}
		for j in (which_idx-1)..(ARROWS - 2*nest_count-2) {
			chooser.out[j] = chooser.out[j+1];
		}

		nest2(nest_count+1, which, chooser, edge_map, vertex_map, genus);

	} else {
		let idx0 = chooser.out[0];
		let idx1 = chooser.out[1];
		edge_map.out[idx0] = idx1;
		edge_map.out[idx1] = idx0;
		count(genus, edge_map, vertex_map);
	}
}

fn nest(nest_count: usize, which: &mut PermutationStructure, chooser: &mut PermutationStructure, 
	edge_map: &mut PermutationStructure, vertex_map: &PermutationStructure, 
	genus: &mut [usize]) {

	if nest_count < ARROWS/2-1 {
		for i in 1..(ARROWS -2*nest_count) {
			which.out[nest_count+1] = i;
			nest(nest_count+1, which, chooser, edge_map, vertex_map, genus);
		}
	} else {
		for j in 0..ARROWS {
			chooser.out[j] = j;
		}
		nest2(0, which, chooser, edge_map, vertex_map, genus);
	}



}

// ----------------------------------------------------------------------

fn main() {
	
	let mut vertex_map = convert_vertex_map();
	let mut edge_map = PermutationStructure {
		out: [0; ARROWS],
	};
	let mut chooser = PermutationStructure {
		out: [0; ARROWS],
	};
	let mut which = PermutationStructure {
		out: [0; ARROWS],
	};
	let mut genus: [usize; MAX_GENUS] = [0; MAX_GENUS];

	nest(0, &mut which, &mut chooser, &mut edge_map, &mut vertex_map, &mut genus);

	println!("Genus Counts:");
	for g in 0..MAX_GENUS {
		println!("Genus {} = {}", g, genus[g] )
	}


}



