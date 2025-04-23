use parry3d_f64::na::Point3;
use tracing::info;

use crate::{
    computational_geometry::intersection_checking,
    error,
    structures::{DFNGen, Poly, PolyOptions},
};

fn in_polygon_boundary(
    num_of_domain_vertices: usize,
    domain_vertices: &[Point3<f64>],
    x: f64,
    y: f64,
) -> bool {
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
    let mut j = num_of_domain_vertices - 1;
    while i < num_of_domain_vertices {
        if ((domain_vertices[i].y > y) != (domain_vertices[j].y > y))
            && (x
                < (domain_vertices[j].x - domain_vertices[i].x) * (y - domain_vertices[i].y)
                    / (domain_vertices[j].y - domain_vertices[i].y)
                    + domain_vertices[i].x)
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

/// Remove Fractures Outside 2D polygon domain
///
/// Function is designed to be used AFTER DFN generation
/// Originally created to compare the difference in distributions
/// if small fractures were removed after a DFN was generated
/// as opposed limiting their insertion during DFN generation.
///
/// The minimum size options in the input files will still be used.
/// If the user wishes to also removeFractures after DFN generation,
/// 'minSize' here must be larger than the minimum size fractures in the
/// input file. Fratrues with radii less than 'minSize' will be removed
/// after DFN has been created.
///
/// # Arguments
///
/// * `h` - Minimum feature size
/// * `eps` - Epsilon value for floating point comparisons
/// * `r_fram` - Uses a relaxed version of the FRAM algorithm. The mesh may not be perfectly conforming
/// * `disable_fram` - If true, FRAM is disabled
/// * `triple_intersections` - If true, triple intersections are accepted
/// * `num_of_domain_vertices` - Number of vertices in the domain polygon
/// * `domain_vertices` -
/// * `accepted_polys` - Array of accepted polygons
/// * `int_pts` - Array of accepted intersections
/// * `triple_points` - Array of all triple intersection points
/// * `pstats` - Stats structure (DFN Statisctics)
#[allow(clippy::too_many_arguments)]
pub fn polygon_boundary(
    opts: &PolyOptions,
    num_of_domain_vertices: usize,
    domain_vertices: &[Point3<f64>],
    dfngen: &mut DFNGen,
) {
    let mut final_poly_list: Vec<Poly> = Vec::new();

    // Clear GroupData
    dfngen.pstats.group_data.clear();
    // Clear FractGroup
    dfngen.pstats.fract_group.clear();
    // Clear Triple Points
    dfngen.triple_points.clear();
    // Clear IntPoints
    dfngen.intpts.clear();
    // Re-init nextGroupNum
    dfngen.pstats.next_group_num = 1;

    for i in 0..dfngen.accepted_poly.len() {
        let x = dfngen.accepted_poly[i].translation[0];
        let y = dfngen.accepted_poly[i].translation[1];

        // cout << "fracture " << i + 1 << " center " << x << "," << y << endl;
        if !in_polygon_boundary(num_of_domain_vertices, domain_vertices, x, y) {
            dfngen.accepted_poly.clear();
            continue;
        }

        let new_poly = &mut dfngen.accepted_poly[i];
        new_poly.group_num = 0; // Reset cluster group number
        new_poly.intersection_index.clear(); // Remove ref to old intersections
                                             // Find line of intersection and FRAM check
        let reject_code = if opts.disable_fram {
            0
        } else {
            intersection_checking(
                opts.h,
                opts.eps,
                opts.r_fram,
                opts.disable_fram,
                opts.triple_intersections,
                new_poly,
                &mut final_poly_list,
                &mut dfngen.intpts,
                &mut dfngen.pstats,
                &mut dfngen.triple_points,
            )
        };

        // IF POLY ACCEPTED:
        if reject_code == 0 {
            // Intersections are ok
            // SAVING POLYGON (intersection and triple points saved witchin intersectionChecking())
            final_poly_list.push(new_poly.clone()); // SAVE newPoly to accepted polys list
        } else {
            // Poly rejected
            error!("Error rebuilding dfn, previously accepted fracture was rejected during DFN rebuild.");
        }
    }

    info!("Rebuilding DFN complete.");
    dfngen.accepted_poly.clear();
    dfngen.accepted_poly.extend(final_poly_list);
}
