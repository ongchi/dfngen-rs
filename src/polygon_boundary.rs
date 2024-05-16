use parry3d::na::Point3;

use crate::{
    computational_geometry::intersection_checking,
    read_input::Input,
    structures::{IntPoints, Poly, Stats},
};

fn in_polygon_boundary(input: &Input, x: f64, y: f64) -> bool {
    // Checks if a point is within the polygon domain using a ray casting algorithm.
    // 1) create a parametric equation for a ray coming from the new Points
    // * r(x) = newPoint.x + t * 1
    // * r(y) = newPoint.y + t * 1
    // * * This ray just goes along the x axis
    // * 2) loop through initial vertices of the polygon.
    // * Create a parametric equation for each edge
    // * l(x) = x + u * mx
    // * l(y) = y + u * my
    // * 3) find point of intersection between the two rays
    // * 4) determine if the point is on the boundary of the polygon.
    // * 5) count the number of times the ray emitting from the point crosses the boundary of the polygon
    // * if the count is zero or even, the point is outside the domain,
    // * if the count is odd, the point is inside the domain
    //
    // Ray Casting adapted from
    //
    // https://wrf.ecse.rpi.edu/Research/Short_Notes/pnpoly.html#The%20Inequality%20Tests%20are%20Tricky
    //
    // License to Use
    //
    // Copyright (c) 1970-2003, Wm. Randolph Franklin
    //
    // Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
    //
    // Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimers.
    // Redistributions in binary form must reproduce the above copyright notice in the documentation and/or other materials provided with the distribution.
    // The name of W. Randolph Franklin may not be used to endorse or promote products derived from this Software without specific prior written permission.
    //
    // THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
    //

    let mut c = false;
    let mut i = 0;
    let mut j = input.numOfDomainVertices - 1;
    while i < input.numOfDomainVertices {
        if ((input.domainVertices[i].y > y) != (input.domainVertices[j].y > y))
            && (x
                < (input.domainVertices[j].x - input.domainVertices[i].x)
                    * (y - input.domainVertices[i].y)
                    / (input.domainVertices[j].y - input.domainVertices[i].y)
                    + input.domainVertices[i].x)
        {
            // flips back and forth between 0 and 1
            // will be 0 for no crossings or an even number of them
            // will be 1 for odd number of crossings
            c = !c;
        }
        j = i;
        i += 1;
    }

    // If the number crossing is odd, then the point is inside of the domain.
    // if the number crossing is zero or even, then the point is outside of the domain.
    c
}

// ***********************************************************************
// ***********  Remove Fractures Outside 2D polygon domain  **************
// Function is designed to be used AFTER DFN generation
// Originally created to compare the difference in distributions
// if small fractures were removed after a DFN was generated
// as opposed limiting their insertion during DFN generation.
//
// The minimum size options in the input files will still be used.
// If the user wishes to also removeFractures after DFN generation,
// 'minSize' here must be larger than the minimum size fractures in the
// input file. Fratrues with radii less than 'minSize' will be removed
// after DFN has been created.
// Arg 1: Minimum size. All fractures smaller than 'minSize' will be
//        removed.
// Arg 2: Array of accepted polygons
// Arg 3: Array of accepted intersections
// Arg 4: Array of all triple intersection points
// Arg 5: Stats structure (DFN Statisctics)
pub fn polygon_boundary(
    input: &Input,
    accepted_polys: &mut Vec<Poly>,
    int_pts: &mut Vec<IntPoints>,
    triple_points: &mut Vec<Point3<f64>>,
    pstats: &mut Stats,
) {
    let mut final_poly_list: Vec<Poly> = Vec::new();

    // Clear GroupData
    pstats.group_data.clear();
    // Clear FractGroup
    pstats.fract_group.clear();
    // Clear Triple Points
    triple_points.clear();
    // Clear IntPoints
    int_pts.clear();
    // Re-init nextGroupNum
    pstats.next_group_num = 1;

    for i in 0..accepted_polys.len() {
        let x = accepted_polys[i].translation[0];
        let y = accepted_polys[i].translation[1];

        // cout << "fracture " << i + 1 << " center " << x << "," << y << endl;
        if !in_polygon_boundary(input, x, y) {
            accepted_polys.clear();
            continue;
        }

        let new_poly = &mut accepted_polys[i];
        new_poly.group_num = 0; // Reset cluster group number
        new_poly.intersection_index.clear(); // Remove ref to old intersections
                                             // Find line of intersection and FRAM check
        let reject_code = if input.disableFram {
            0
        } else {
            intersection_checking(
                input,
                new_poly,
                &mut final_poly_list,
                int_pts,
                pstats,
                triple_points,
            )
        };

        // IF POLY ACCEPTED:
        if reject_code == 0 {
            // Intersections are ok
            // SAVING POLYGON (intersection and triple points saved witchin intersectionChecking())
            final_poly_list.push(new_poly.clone()); // SAVE newPoly to accepted polys list
        } else {
            // Poly rejected
            println!("Error rebuilding dfn, previously accepted fracture was rejected during DFN rebuild.");
        }
    }

    println!("Rebuilding DFN complete.");
    accepted_polys.clear();
    accepted_polys.extend(final_poly_list);
}
