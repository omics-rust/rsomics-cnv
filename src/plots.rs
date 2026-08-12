use std::fmt;
use std::ops::Range;

use plotters::coord::Shift;
use plotters::prelude::*;
use rsomics_common::{Result, RsomicsError};

use crate::call::{CallResult, ChromosomeCall};
use crate::polysomy::{ChromosomePolysomy, PolysomyResult};

const MAX_POINTS: usize = 10_000;

pub(crate) fn for_each_call_plot(
    result: &CallResult,
    minimum_quality: f64,
    mut output: impl FnMut(String, String) -> Result<()>,
) -> Result<()> {
    for chromosome in &result.chromosomes {
        if let Some((name, svg)) = call_plot(
            chromosome,
            &result.sample,
            result.control_sample.as_deref(),
            minimum_quality,
        )? {
            output(name, svg)?;
        }
    }
    Ok(())
}

pub(crate) fn call_plot(
    chromosome: &ChromosomeCall,
    query: &str,
    control: Option<&str>,
    minimum_quality: f64,
) -> Result<Option<(String, String)>> {
    if !minimum_quality.is_finite() || minimum_quality < 0.0 {
        return Err(invalid(
            "plot threshold must be finite and greater than or equal to zero",
        ));
    }
    let maximum = chromosome
        .regions
        .iter()
        .map(|region| region.quality)
        .fold(f64::NEG_INFINITY, f64::max);
    if maximum < minimum_quality {
        return Ok(None);
    }
    let name = match control {
        Some(control) => format!(
            "plot.{}.{}.{}.svg",
            filename_component(control),
            filename_component(query),
            filename_component(&chromosome.reference_name)
        ),
        None => format!(
            "plot.{}.{}.svg",
            filename_component(query),
            filename_component(&chromosome.reference_name)
        ),
    };
    Ok(Some((name, render_call(chromosome, query, control)?)))
}

pub(crate) fn for_each_polysomy_plot(
    result: &PolysomyResult,
    mut output: impl FnMut(String, String) -> Result<()>,
) -> Result<()> {
    for chromosome in &result.chromosomes {
        output(
            format!(
                "distribution.{}.svg",
                filename_component(&chromosome.distribution.reference_name)
            ),
            render_distribution(chromosome)?,
        )?;
    }
    output("copy-number.svg".to_owned(), render_copy_numbers(result)?)
}

fn render_call(chromosome: &ChromosomeCall, query: &str, control: Option<&str>) -> Result<String> {
    let samples = if let Some(control) = control {
        vec![(control, true, RED), (query, false, BLUE)]
    } else {
        vec![(query, false, BLUE)]
    };
    let panels = samples.len() * 3;
    let mut svg = String::new();
    {
        let root =
            SVGBackend::with_string(&mut svg, (960, panels as u32 * 230)).into_drawing_area();
        root.fill(&WHITE).map_err(plot_error)?;
        let areas = root.split_evenly((panels, 1));
        let x = coordinate_range(chromosome)?;
        for (sample_index, (sample, control, color)) in samples.into_iter().enumerate() {
            let measurements = sampled_measurements(chromosome, control);
            let lrr = measurements
                .iter()
                .filter_map(|(position, measurement)| {
                    measurement.lrr.map(|value| (*position, value))
                })
                .collect::<Vec<_>>();
            let baf = measurements
                .iter()
                .filter_map(|(position, measurement)| {
                    measurement.baf.map(|value| (*position, value))
                })
                .collect::<Vec<_>>();
            scatter_panel(
                &areas[sample_index * 3],
                &format!("{sample} LRR"),
                "LRR",
                x.clone(),
                value_range(&lrr, -1.0..1.0),
                &lrr,
                color,
            )?;
            scatter_panel(
                &areas[sample_index * 3 + 1],
                &format!("{sample} BAF"),
                "BAF",
                x.clone(),
                -0.05..1.05,
                &baf,
                color,
            )?;
            copy_number_panel(
                &areas[sample_index * 3 + 2],
                &format!("{sample} Copy number — {}", chromosome.reference_name),
                x.clone(),
                chromosome,
                control,
                color,
            )?;
        }
        root.present().map_err(plot_error)?;
    }
    Ok(svg)
}

fn scatter_panel<DB: DrawingBackend>(
    area: &DrawingArea<DB, Shift>,
    caption: &str,
    y_label: &str,
    x: Range<u32>,
    y: Range<f64>,
    points: &[(u32, f64)],
    color: RGBColor,
) -> Result<()>
where
    DB::ErrorType: fmt::Debug,
{
    let mut chart = ChartBuilder::on(area)
        .caption(caption, ("sans-serif", 18))
        .margin(6)
        .x_label_area_size(32)
        .y_label_area_size(52)
        .build_cartesian_2d(x, y)
        .map_err(plot_error)?;
    chart
        .configure_mesh()
        .x_desc("Position")
        .y_desc(y_label)
        .light_line_style(WHITE)
        .draw()
        .map_err(plot_error)?;
    chart
        .draw_series(
            points
                .iter()
                .map(|&(position, value)| Circle::new((position, value), 2, color.filled())),
        )
        .map_err(plot_error)?;
    Ok(())
}

fn copy_number_panel<DB: DrawingBackend>(
    area: &DrawingArea<DB, Shift>,
    caption: &str,
    x: Range<u32>,
    chromosome: &ChromosomeCall,
    control: bool,
    color: RGBColor,
) -> Result<()>
where
    DB::ErrorType: fmt::Debug,
{
    let mut chart = ChartBuilder::on(area)
        .caption(caption, ("sans-serif", 18))
        .margin(6)
        .x_label_area_size(32)
        .y_label_area_size(52)
        .build_cartesian_2d(x, -0.25f64..4.25f64)
        .map_err(plot_error)?;
    chart
        .configure_mesh()
        .x_desc("Position")
        .y_desc("Copy number")
        .y_labels(5)
        .light_line_style(WHITE)
        .draw()
        .map_err(plot_error)?;
    let points = copy_number_steps(chromosome, control)?;
    chart
        .draw_series(LineSeries::new(
            points.iter().copied(),
            color.stroke_width(2),
        ))
        .map_err(plot_error)?;
    chart
        .draw_series(
            points
                .iter()
                .map(|&(position, value)| Circle::new((position, value), 2, color.filled())),
        )
        .map_err(plot_error)?;
    Ok(())
}

fn render_distribution(chromosome: &ChromosomePolysomy) -> Result<String> {
    let distribution = &chromosome.distribution;
    let decision = if chromosome.copy_number >= 0.0 {
        format!("CN {:.2}", chromosome.copy_number)
    } else {
        "unresolved".to_owned()
    };
    let selected = chromosome
        .candidates
        .iter()
        .find(|candidate| candidate.selected);
    let maximum = distribution
        .bins
        .iter()
        .map(|bin| bin.normalized_count)
        .chain(
            selected
                .into_iter()
                .flat_map(|candidate| &candidate.curves)
                .flat_map(|curve| &curve.fitted)
                .map(|bin| bin.normalized_count),
        )
        .fold(1.0f64, f64::max)
        * 1.08;
    let mut svg = String::new();
    {
        let root = SVGBackend::with_string(&mut svg, (900, 620)).into_drawing_area();
        root.fill(&WHITE).map_err(plot_error)?;
        let mut chart = ChartBuilder::on(&root)
            .caption(
                format!(
                    "BAF distribution — {} — {decision}",
                    distribution.reference_name
                ),
                ("sans-serif", 22),
            )
            .margin(15)
            .x_label_area_size(45)
            .y_label_area_size(55)
            .build_cartesian_2d(0.0f64..1.0f64, 0.0f64..maximum)
            .map_err(plot_error)?;
        chart
            .configure_mesh()
            .x_desc("BAF")
            .y_desc("Normalized frequency")
            .light_line_style(WHITE)
            .draw()
            .map_err(plot_error)?;
        chart
            .draw_series(LineSeries::new(
                distribution
                    .bins
                    .iter()
                    .map(|bin| (bin.baf, bin.normalized_count)),
                BLACK.stroke_width(2),
            ))
            .map_err(plot_error)?
            .label("Distribution")
            .legend(|(x, y)| PathElement::new([(x, y), (x + 24, y)], BLACK.stroke_width(2)));
        if let Some(candidate) = selected {
            for (index, curve) in candidate.curves.iter().enumerate() {
                let annotation = chart
                    .draw_series(LineSeries::new(
                        curve
                            .fitted
                            .iter()
                            .map(|bin| (bin.baf, bin.normalized_count)),
                        GREEN.stroke_width(2),
                    ))
                    .map_err(plot_error)?;
                if index == 0 {
                    annotation.label("Selected fit").legend(|(x, y)| {
                        PathElement::new([(x, y), (x + 24, y)], GREEN.stroke_width(2))
                    });
                }
            }
        }
        chart
            .configure_series_labels()
            .background_style(WHITE.mix(0.9))
            .border_style(BLACK)
            .draw()
            .map_err(plot_error)?;
        root.present().map_err(plot_error)?;
    }
    Ok(svg)
}

fn render_copy_numbers(result: &PolysomyResult) -> Result<String> {
    let mut svg = String::new();
    {
        let root = SVGBackend::with_string(&mut svg, (900, 560)).into_drawing_area();
        root.fill(&WHITE).map_err(plot_error)?;
        let count = result.chromosomes.len();
        let names = result
            .chromosomes
            .iter()
            .map(|chromosome| chromosome.distribution.reference_name.as_str())
            .collect::<Vec<_>>();
        let mut chart = ChartBuilder::on(&root)
            .caption(
                format!("Copy number — {}", result.sample),
                ("sans-serif", 22),
            )
            .margin(15)
            .x_label_area_size(55)
            .y_label_area_size(55)
            .build_cartesian_2d(-0.5f64..count as f64 - 0.5, -0.25f64..5.0f64)
            .map_err(plot_error)?;
        chart
            .configure_mesh()
            .x_desc("Chromosome")
            .y_desc("Copy number")
            .x_labels(count)
            .x_label_formatter(&|value| {
                let index = value.round() as isize;
                usize::try_from(index)
                    .ok()
                    .and_then(|index| names.get(index))
                    .copied()
                    .unwrap_or("")
                    .to_owned()
            })
            .light_line_style(WHITE)
            .draw()
            .map_err(plot_error)?;
        chart
            .draw_series(
                result
                    .chromosomes
                    .iter()
                    .enumerate()
                    .filter(|(_, chromosome)| chromosome.copy_number >= 0.0)
                    .map(|(index, chromosome)| {
                        Circle::new((index as f64, chromosome.copy_number), 5, RED.filled())
                    }),
            )
            .map_err(plot_error)?;
        let unresolved = result
            .chromosomes
            .iter()
            .enumerate()
            .filter(|(_, chromosome)| chromosome.copy_number < 0.0)
            .map(|(index, _)| (index as f64, 0.0))
            .collect::<Vec<_>>();
        if !unresolved.is_empty() {
            chart
                .draw_series(
                    unresolved
                        .into_iter()
                        .map(|point| Cross::new(point, 6, BLACK.stroke_width(2))),
                )
                .map_err(plot_error)?
                .label("Unresolved")
                .legend(|(x, y)| Cross::new((x + 12, y), 5, BLACK.stroke_width(2)));
            chart
                .configure_series_labels()
                .background_style(WHITE.mix(0.9))
                .border_style(BLACK)
                .draw()
                .map_err(plot_error)?;
        }
        root.present().map_err(plot_error)?;
    }
    Ok(svg)
}

fn coordinate_range(chromosome: &ChromosomeCall) -> Result<Range<u32>> {
    let start = chromosome
        .sites
        .first()
        .ok_or_else(|| invalid("cannot plot a chromosome without sites"))?
        .position;
    let end = chromosome
        .sites
        .last()
        .ok_or_else(|| invalid("cannot plot a chromosome without sites"))?
        .position;
    Ok(if start == end {
        start.saturating_sub(1)..end.saturating_add(1).max(end)
    } else {
        start..end
    })
}

fn sampled_measurements(
    chromosome: &ChromosomeCall,
    control: bool,
) -> Vec<(u32, crate::signals::Measurement)> {
    sampled_indices(chromosome.sites.len())
        .into_iter()
        .filter_map(|index| {
            let site = &chromosome.sites[index];
            let measurement = if control {
                site.control_measurement?
            } else {
                site.measurement
            };
            Some((site.position, measurement))
        })
        .collect()
}

fn sampled_indices(length: usize) -> Vec<usize> {
    if length <= MAX_POINTS {
        return (0..length).collect();
    }
    (0..MAX_POINTS)
        .map(|index| index * (length - 1) / (MAX_POINTS - 1))
        .collect()
}

fn value_range(points: &[(u32, f64)], fallback: Range<f64>) -> Range<f64> {
    let minimum = points
        .iter()
        .map(|(_, value)| *value)
        .fold(f64::INFINITY, f64::min);
    let maximum = points
        .iter()
        .map(|(_, value)| *value)
        .fold(f64::NEG_INFINITY, f64::max);
    if !minimum.is_finite() || !maximum.is_finite() {
        return fallback;
    }
    let span = (maximum - minimum).max(0.1);
    minimum - span * 0.08..maximum + span * 0.08
}

fn copy_number_steps(chromosome: &ChromosomeCall, control: bool) -> Result<Vec<(u32, f64)>> {
    let mut points = Vec::with_capacity(chromosome.regions.len() * 3);
    for (index, region) in chromosome.regions.iter().enumerate() {
        let copy_number = if control {
            region
                .control_copy_number
                .ok_or_else(|| invalid("control plot is missing a copy-number state"))?
        } else {
            region.copy_number
        };
        let value = f64::from(copy_number);
        points.push((region.start, value));
        points.push((region.end, value));
        if let Some(next) = chromosome.regions.get(index + 1) {
            points.push((next.start, value));
        }
    }
    Ok(points)
}

fn filename_component(value: &str) -> String {
    let mut output = String::new();
    for byte in value.bytes() {
        if byte.is_ascii_alphanumeric() || matches!(byte, b'.' | b'-' | b'_') {
            output.push(char::from(byte));
        } else {
            output.push_str(&format!("%{byte:02X}"));
        }
    }
    output
}

fn plot_error(error: impl fmt::Debug) -> RsomicsError {
    invalid(format!("rendering SVG plot: {error:?}"))
}

fn invalid(message: impl Into<String>) -> RsomicsError {
    RsomicsError::InvalidInput(message.into())
}
