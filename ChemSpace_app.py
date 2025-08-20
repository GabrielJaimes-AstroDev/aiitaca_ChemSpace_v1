import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import json
import re
import os
from io import StringIO
import tempfile

# Page configuration
st.set_page_config(
    page_title="Chemical Space Visualization",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        color: #1f77b4;
        text-align: center;
        margin-bottom: 2rem;
    }
    .sub-header {
        font-size: 1.5rem;
        color: #2ca02c;
        margin-bottom: 1rem;
    }
    .info-box {
        background-color: #f0f8ff;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #1f77b4;
        margin-bottom: 1rem;
    }
    .error-box {
        background-color: #ffebee;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #d32f2f;
        margin-bottom: 1rem;
    }
    .success-box {
        background-color: #e8f5e8;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #2ca02c;
        margin-bottom: 1rem;
    }
</style>
""", unsafe_allow_html=True)

def repair_json_file(input_path, output_path):
    """Repair a corrupted JSON file with multiple error handling strategies"""
    try:
        with open(input_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        st.info(f"Original file size: {len(content)} characters")
        
        # Strategy 1: Try direct JSON loading first
        try:
            data = json.loads(content)
            with open(output_path, 'w', encoding='utf-8') as f:
                json.dump(data, f, indent=2)
            return True, "JSON loaded successfully without repair"
        except json.JSONDecodeError as e:
            st.warning(f"Initial JSON parse failed: {e}")
        
        # Strategy 2: Remove problematic lines and try to salvage data
        lines = content.split('\n')
        st.info(f"Total lines: {len(lines)}")
        
        # Find and remove problematic lines
        repaired_lines = []
        removed_lines = 0
        
        for i, line in enumerate(lines, 1):
            line = line.strip()
            if not line:
                continue
                
            # Skip lines that are clearly malformed
            if line.endswith(',') and not line.endswith('],') and not line.endswith('},'):
                st.warning(f"Removing malformed line {i}: {line[:100]}...")
                removed_lines += 1
                continue
                
            # Fix common array issues
            if line.endswith(','):
                line = line.rstrip(',')
            
            repaired_lines.append(line)
        
        repaired_content = '\n'.join(repaired_lines)
        st.info(f"Removed {removed_lines} problematic lines")
        
        # Strategy 3: Extract data using regex patterns
        def extract_array(pattern, content):
            matches = re.findall(pattern, content, re.DOTALL)
            if matches:
                return matches[0]
            return None
        
        # Extract points array
        points_pattern = r'"points"\s*:\s*\[(.*?)\]'
        points_content = extract_array(points_pattern, content)
        
        # Extract other arrays
        databases_pattern = r'"databases"\s*:\s*\[(.*?)\]'
        labels_pattern = r'"labels"\s*:\s*\[(.*?)\]'
        formulas_pattern = r'"formulas"\s*:\s*\[(.*?)\]'
        colors_pattern = r'"colors"\s*:\s*\[(.*?)\]'
        
        points_content = extract_array(points_pattern, content)
        databases_content = extract_array(databases_pattern, content)
        labels_content = extract_array(labels_pattern, content)
        formulas_content = extract_array(formulas_pattern, content)
        colors_content = extract_array(colors_pattern, content)
        
        if all([points_content, databases_content, labels_content, formulas_content, colors_content]):
            # Parse arrays individually
            def parse_array(array_str, is_string=True):
                try:
                    if is_string:
                        # Extract string values
                        values = re.findall(r'"([^"]*)"', array_str)
                        return values
                    else:
                        # Extract numeric arrays
                        arrays = re.findall(r'\[([^\]]*)\]', array_str)
                        result = []
                        for arr in arrays:
                            numbers = re.findall(r'[-+]?\d*\.\d+|\d+', arr)
                            result.append([float(x) for x in numbers])
                        return result
                except:
                    return []
            
            points = parse_array(points_content, False)
            databases = parse_array(databases_content, True)
            labels = parse_array(labels_content, True)
            formulas = parse_array(formulas_content, True)
            colors = parse_array(colors_content, True)
            
            # Ensure all arrays have the same length
            min_length = min(len(points), len(databases), len(labels), len(formulas), len(colors))
            
            data = {
                'points': points[:min_length],
                'databases': databases[:min_length],
                'labels': labels[:min_length],
                'formulas': formulas[:min_length],
                'colors': colors[:min_length]
            }
            
            with open(output_path, 'w', encoding='utf-8') as f:
                json.dump(data, f, indent=2)
            
            return True, f"Successfully extracted {min_length} molecules using regex parsing"
        
        # Strategy 4: Try to parse as NDJSON (newline-delimited JSON)
        try:
            data_lines = []
            for line in repaired_lines:
                if line.strip() and not line.strip().startswith('//'):
                    try:
                        data = json.loads(line)
                        data_lines.append(data)
                    except:
                        continue
            
            if data_lines:
                # If we have NDJSON data, try to combine it
                combined_data = {}
                for line_data in data_lines:
                    for key, value in line_data.items():
                        if key not in combined_data:
                            combined_data[key] = []
                        if isinstance(value, list):
                            combined_data[key].extend(value)
                
                with open(output_path, 'w', encoding='utf-8') as f:
                    json.dump(combined_data, f, indent=2)
                
                return True, "Successfully parsed as NDJSON format"
                
        except Exception as e:
            st.warning(f"NDJSON parsing failed: {e}")
        
        return False, "All repair strategies failed"
        
    except Exception as e:
        return False, f"File processing error: {str(e)}"

def create_sample_data():
    """Create sample data for demonstration"""
    return {
        'points': [[0.1, 0.2], [0.3, 0.4], [0.5, 0.6]],
        'databases': ['GUAPOS', 'PubChem', 'All_Discoveries'],
        'labels': ['Molecule1', 'Molecule2', 'Molecule3'],
        'formulas': ['H2O', 'C6H12O6', 'NaCl'],
        'colors': ['#1f77b4', '#9467bd', '#2ca02c'],
        'smiles': ['O', 'C(C1C(C(C(C1O)O)O)O)O', '[Na+].[Cl-]'],
        'detected_status': [1, 0, 1],
        'knn_connections': [
            {'guapos_index': 0, 'neighbor_index': 1, 'distance': 0.15, 'rank': 1}
        ]
    }

def format_chemical_formula(formula):
    """Format chemical formula with subscripts"""
    if not formula or pd.isna(formula) or formula == 'Unknown':
        return ""
    
    formula = re.sub(r'(\d+)', r'<sub>\1</sub>', str(formula))
    formula = re.sub(r'([+-]\d*)', r'<sup>\1</sup>', formula)
    
    return formula

def create_interactive_plot(plot_data):
    """Create interactive Plotly visualization"""
    try:
        # Create DataFrame from plot data
        df_points = pd.DataFrame({
            'x': [point[0] for point in plot_data['points']],
            'y': [point[1] for point in plot_data['points']],
            'label': plot_data['labels'],
            'database': plot_data['databases'],
            'formula': plot_data['formulas'],
            'color': plot_data['colors'],
            'index': range(len(plot_data['points']))
        })
        
        fig = go.Figure()
        
        # Add points for each database
        databases = df_points['database'].unique()
        for db in databases:
            db_data = df_points[df_points['database'] == db]
            
            hover_text = []
            for _, row in db_data.iterrows():
                text = f"Database: {row['database']}<br>"
                if row['label']:
                    text += f"Name: {row['label']}<br>"
                if row['formula']:
                    formatted_formula = format_chemical_formula(row['formula'])
                    text += f"Formula: {formatted_formula}"
                hover_text.append(text)
            
            fig.add_trace(go.Scatter(
                x=db_data['x'],
                y=db_data['y'],
                mode='markers',
                marker=dict(
                    color=db_data['color'].iloc[0],
                    size=12,
                    opacity=0.7
                ),
                name=db,
                hovertext=hover_text,
                hoverinfo='text'
            ))
        
        fig.update_layout(
            title="Chemical Space Visualization",
            xaxis_title="Dimension 1",
            yaxis_title="Dimension 2",
            hovermode='closest',
            height=600
        )
        
        fig.update_xaxes(showgrid=False, zeroline=False, showticklabels=False)
        fig.update_yaxes(showgrid=False, zeroline=False, showticklabels=False)
        
        return fig
        
    except Exception as e:
        st.error(f"Error creating plot: {str(e)}")
        return go.Figure()

def main():
    st.markdown('<h1 class="main-header">🧪 Chemical Space Visualization</h1>', unsafe_allow_html=True)
    
    # File upload section
    st.markdown("---")
    st.header("📁 Upload JSON File")
    
    uploaded_file = st.file_uploader(
        "Upload your plot_data.json file",
        type=['json'],
        help="Select the JSON file generated by the analysis script"
    )
    
    use_sample = st.checkbox("Use sample data instead", value=False)
    
    plot_data = None
    
    if use_sample:
        plot_data = create_sample_data()
        st.success("✅ Using sample data for demonstration!")
    
    elif uploaded_file is not None:
        try:
            # Save uploaded file
            with tempfile.NamedTemporaryFile(delete=False, suffix='.json') as tmp_file:
                tmp_file.write(uploaded_file.getvalue())
                input_path = tmp_file.name
            
            # Create output file
            output_path = "repaired_plot_data.json"
            
            # Attempt to repair the JSON
            with st.spinner("🛠️ Repairing JSON file..."):
                success, message = repair_json_file(input_path, output_path)
            
            if success:
                st.success(f"✅ {message}")
                
                # Load the repaired data
                try:
                    with open(output_path, 'r', encoding='utf-8') as f:
                        plot_data = json.load(f)
                    st.success("✅ Repaired data loaded successfully!")
                    
                    # Show data summary
                    st.info(f"""
                    **Data Summary:**
                    - Molecules: {len(plot_data['points'])}
                    - Databases: {len(set(plot_data['databases']))}
                    - Labels: {len(plot_data['labels'])}
                    - Formulas: {len(plot_data['formulas'])}
                    """)
                    
                except Exception as e:
                    st.error(f"❌ Error loading repaired file: {str(e)}")
                    
            else:
                st.error(f"❌ {message}")
                st.info("""
                **Unable to repair the JSON file. Possible reasons:**
                1. The file is severely corrupted
                2. The file format is not as expected
                3. The file is too large or complex
                
                **Please try:**
                - Generating a new JSON file from your analysis script
                - Checking the file format matches the expected structure
                - Using the sample data option to verify the app works
                """)
                
        except Exception as e:
            st.error(f"❌ File processing error: {str(e)}")
    
    else:
        st.info("""
        👆 **Please upload a JSON file or use the sample data option**
        
        **Expected JSON structure:**
        ```json
        {
          "points": [[x1, y1], [x2, y2], ...],
          "databases": ["db1", "db2", ...],
          "labels": ["name1", "name2", ...],
          "formulas": ["formula1", "formula2", ...],
          "colors": ["color1", "color2", ...]
        }
        ```
        """)
    
    # Display visualization if we have data
    if plot_data:
        st.markdown("---")
        st.header("📊 Visualization")
        
        # Create and display plot
        fig = create_interactive_plot(plot_data)
        st.plotly_chart(fig, use_container_width=True)
        
        # Data statistics
        st.markdown("---")
        st.header("📈 Statistics")
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.metric("Total Molecules", len(plot_data['points']))
        
        with col2:
            unique_dbs = len(set(plot_data['databases']))
            st.metric("Unique Databases", unique_dbs)
        
        with col3:
            guapos_count = plot_data['databases'].count('GUAPOS') if 'GUAPOS' in plot_data['databases'] else 0
            st.metric("GUAPOS Molecules", guapos_count)
        
        # Database distribution
        db_counts = pd.Series(plot_data['databases']).value_counts()
        st.write("**Database Distribution:**")
        for db, count in db_counts.items():
            st.write(f"- {db}: {count} molecules")
    
    # Add debug information
    st.markdown("---")
    st.header("🛠️ Debug Tools")
    
    if uploaded_file and not use_sample:
        st.write("**Uploaded File Information:**")
        st.write(f"- Name: {uploaded_file.name}")
        st.write(f"- Size: {uploaded_file.size} bytes")
        st.write(f"- Type: {uploaded_file.type}")
        
        if st.button("🔄 Try Alternative Repair Method"):
            with st.spinner("Trying alternative repair method..."):
                try:
                    # Try a different approach: read line by line
                    content = uploaded_file.getvalue().decode('utf-8')
                    lines = content.split('\n')
                    
                    st.info(f"File has {len(lines)} lines")
                    st.info(f"First few lines:")
                    for i, line in enumerate(lines[:5]):
                        st.write(f"Line {i+1}: {line[:100]}...")
                    
                    # Look for specific patterns
                    points_lines = [i for i, line in enumerate(lines) if '"points"' in line]
                    if points_lines:
                        st.success(f"Found 'points' array starting at line {points_lines[0] + 1}")
                    
                except Exception as e:
                    st.error(f"Alternative repair failed: {e}")

if __name__ == "__main__":
    main()
