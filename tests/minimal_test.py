import os
import sys

def add_module_path():
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    build_dir = os.path.join(project_root, 'build')
    
    if os.path.exists(build_dir):
        print(f"Adding build directory to path: {build_dir}")
        sys.path.insert(0, build_dir)
    else:
        print(f"Build directory not found: {build_dir}")
        return False
    return True

def test_import():
    try:
        from dr_sasa_py import DrSASA, ComputeBackend
        print("Successfully imported dr_sasa_py")
        
        # Try creating an instance
        calculator = DrSASA()
        print("Successfully created DrSASA instance")
        
        return True
    except Exception as e:
        print(f"Error: {e}")
        print(f"Error type: {type(e)}")
        print(f"Error args: {e.args}")
        return False

if __name__ == '__main__':
    if add_module_path():
        if test_import():
            print("All tests passed!")
        else:
            print("Tests failed!")
            sys.exit(1)
    else:
        print("Failed to set up module path!")
        sys.exit(1)