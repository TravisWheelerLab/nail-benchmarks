pub(crate) fn label(index: Option<&str>, name: Option<&str>) -> String {
    match (index, name) {
        (Some(index), Some(name)) => format!("[{index}]({name})"),
        (Some(index), None) => format!("[{index}]"),
        (None, Some(name)) => name.to_string(),
        (None, None) => String::new(),
    }
}

pub(crate) fn filename(index: usize, name: Option<&str>) -> String {
    match name {
        Some(name) => format!("{index}-{}", safe(name)),
        None => index.to_string(),
    }
}

fn safe(name: &str) -> String {
    name.chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || "._-".contains(c) {
                c
            } else {
                '_'
            }
        })
        .collect()
}
