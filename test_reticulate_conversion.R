# Test what reticulate passes to Python for different R objects

library(reticulate)

# Configure Python
use_python("C:/Users/oluwa/OneDrive - Texas Tech University/GlypPRM/Package/test_env/Scripts/python.exe",
           required = TRUE)

# Create test Python function
py_run_string("
def test_fragment_types(fragment_types):
    print(f'Type: {type(fragment_types)}')
    print(f'Value: {fragment_types!r}')
    print(f'Length: {len(fragment_types)}')
    if fragment_types:
        print(f'Is truthy: True')
        print(f'First item: {fragment_types[0]!r}')
        if len(fragment_types) > 0:
            print(f'Iteration:')
            for i, item in enumerate(fragment_types):
                print(f'  [{i}] = {item!r}')
    else:
        print(f'Is falsy!')
    return None
")

test_fn <- py$test_fragment_types

cat("\n=== Test 1: c('BY') directly ===\n")
test_fn(c("BY"))

cat("\n=== Test 2: as.list(c('BY')) ===\n")
test_fn(as.list(c("BY")))

cat("\n=== Test 3: list('BY') ===\n")
test_fn(list("BY"))

cat("\n=== Test 4: c('BY', 'CZ') ===\n")
test_fn(c("BY", "CZ"))

cat("\n=== Test 5: as.list(c('BY', 'CZ')) ===\n")
test_fn(as.list(c("BY", "CZ")))
