#![forbid(unsafe_code)]

mod allele_frequency;
pub mod call;
mod cli;
pub mod emission;
mod fitting;
pub mod hmm;
mod plots;
pub mod polysomy;
pub mod reports;
pub mod selection;
pub mod signals;

#[must_use]
pub fn run_binary() -> std::process::ExitCode {
    cli::run()
}
