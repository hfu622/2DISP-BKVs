//! Export original f64 contours and XML-legal transformations for independent QA.
use crate::{geometry::Geometry, input::LoadedInput, model::*};
use anyhow::Result;
use serde_json::{Value, json};
pub fn output(input: &LoadedInput, g: &Geometry, layout: &Layout) -> Result<Value> {
    g.validate(layout)?;
    output_for_independent_verification(input,g,layout)
}
/// MIP feasibility tolerances can exceed the accelerated oracle's strict contact
/// gate. Export original contours without snapping, so the independent Python
/// verifier can quantify and accept/reject the user's declared tolerances.
pub fn output_for_independent_verification(input:&LoadedInput,g:&Geometry,layout:&Layout)->Result<Value>{
    let mut embedded = serde_json::to_value(&input.instance)?;
    for item in embedded["items"].as_array_mut().unwrap() {
        let id = item["id"].as_u64().unwrap();
        item["shape"]["data"] = serde_json::to_value(&input.exact_vertices[&id])?;
    }
    let placed: Vec<_> = layout
        .placements
        .iter()
        .map(|p| {
            let o = g.orientation(p);
            json!({"item_id":g.types[g.pieces[p.piece]].id,
            "transformation":{"rotation":o.angle.to_radians(),
            "translation":[p.x-o.min_x,p.y-o.min_y]}})
        })
        .collect();
    Ok(
        json!({"instance":embedded,"reproduction_configuration":crate::config::record(),"solution":{"strip_width":layout.length,
        "density":g.total_area/(g.height*layout.length),"layout":{"placed_items":placed}}}),
    )
}
