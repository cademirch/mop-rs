use assert_cmd::prelude::*; 
use std::fs;
use std::io::Read;
use std::process::Command;
use std::path::Path;

#[test]
fn test_moprs_output() -> Result<(), Box<dyn std::error::Error>> {
    
    let mut cmd = Command::cargo_bin("moprs")?;
    
    
    cmd.arg("--d4").arg("tests/data/test.per-base.d4")
        .arg("-m").arg("3")
        .arg("-M").arg("20")
        .arg("-d").arg("1.0")
        .arg("-u").arg("9");

    
    let output_file = "test.bed";

    
    let output = cmd.output()?;
    
    
    fs::write(output_file, &output.stdout)?;

    
    let mut generated_output = String::new();
    fs::File::open(output_file)?
        .read_to_string(&mut generated_output)?;

    
    let expected_output_path = Path::new("tests/data/truth.bed");
    let mut expected_output = String::new();
    fs::File::open(expected_output_path)?
        .read_to_string(&mut expected_output)?;

    
    assert_eq!(generated_output, expected_output);

    // Clean up by removing the generated output file
    fs::remove_file(output_file)?;

    Ok(())
}
