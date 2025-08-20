import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import json
import re
import os
from io import StringIO

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
    .molecule-card {
        background-color: #f8f9fa;
        padding: 0.5rem;
        border-radius: 0.3rem;
        margin: 0.2rem 0;
        border-left: 3px solid #9467bd;
    }
    .neighbor-point {
        background-color: #ffebee;
        padding: 0.5rem;
        border-radius: 0.3rem;
        margin: 0.2rem 0;
        border-left: 3px solid #d32f2f;
    }
    .guapos-detected {
        background-color: #e3f2fd;
        padding: 0.5rem;
        border-radius: 0.3rem;
        margin: 0.2rem 0;
        border-left: 3px solid #1976d2;
    }
    .guapos-not-detected {
        background-color: #f5f5f5;
        padding: 0.5rem;
        border-radius: 0.3rem;
        margin: 0.2rem 0;
        border-left: 3px solid #9e9e9e;
    }
    .error-box {
        background-color: #ffebee;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #d32f2f;
        margin-bottom: 1rem;
    }
</style>
""", unsafe_allow_html=True)

def validate_and_repair_json(file_path):
    """Validate and repair JSON file with detailed error reporting"""
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # Try to parse the JSON
        try:
            data = json.loads(content)
            return data, True, "JSON is valid"
        except json.JSONDecodeError as e:
            st.warning(f"JSON validation error: {e}")
            
            # Try to fix common JSON issues
            try:
                # Remove trailing commas
                content = re.sub(r',\s*}', '}', content)
                content = re.sub(r',\s*]', ']', content)
                
                # Fix unquoted keys
                content = re.sub(r'(\w+):', r'"\1":', content)
                
                # Try to parse again
                data = json.loads(content)
                return data, True, "JSON repaired successfully"
                
            except json.JSONDecodeError as e2:
                # If repair fails, try to extract what we can
                st.error(f"JSON repair failed: {e2}")
                
                # Try to find the approximate location of the error
                lines = content.split('\n')
                if e.lineno - 1 < len(lines):
                    error_line = lines[e.lineno - 1]
                    st.error(f"Error around line {e.lineno}: {error_line}")
                
                return None, False, f"JSON repair failed: {e2}"
                
    except Exception as e:
        return None, False, f"File reading error: {e}"

def load_plot_data(file_path):
    """Load plot data from JSON file with robust error handling"""
    try:
        if not os.path.exists(file_path):
            return None, "File not found"
        
        file_size = os.path.getsize(file_path)
        st.info(f"File size: {file_size / (1024*1024):.2f} MB")
        
        # Validate and repair JSON
        data, is_valid, message = validate_and_repair_json(file_path)
        
        if is_valid:
            # Validate required fields
            required_fields = ['points', 'databases', 'labels', 'formulas', 'colors']
            missing_fields = [field for field in required_fields if field not in data]
            
            if missing_fields:
                return None, f"Missing required fields: {missing_fields}"
            
            # Validate array lengths
            array_lengths = {
                'points': len(data['points']),
                'databases': len(data['databases']),
                'labels': len(data['labels']),
                'formulas': len(data['formulas']),
                'colors': len(data['colors'])
            }
            
            # Check if all arrays have the same length
            lengths = list(array_lengths.values())
            if len(set(lengths)) > 1:
                st.warning(f"Array length mismatch: {array_lengths}")
                # Use the minimum length to avoid index errors
                min_length = min(lengths)
                for key in required_fields:
                    if len(data[key]) > min_length:
                        data[key] = data[key][:min_length]
            
            return data, "Data loaded successfully"
        else:
            return None, message
            
    except Exception as e:
        return None, f"Error loading plot data: {str(e)}"

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
    
    # Handle subscripts
    formula = re.sub(r'(\d+)', r'<sub>\1</sub>', str(formula))
    
    # Handle ionic charges
    formula = re.sub(r'([+-]\d*)', r'<sup>\1</sup>', formula)
    
    return formula

def create_interactive_plot(plot_data):
    """Create interactive Plotly visualization with neighbor highlighting"""
    
    # Create DataFrame from plot data
    df_points = pd.DataFrame({
        'x': [point[0] for point in plot_data['points']],
        'y': [point[1] for point in plot_data['points']],
        'label': plot_data['labels'],
        'database': plot_data['databases'],
        'formula': plot_data['formulas'],
        'color': plot_data['colors'],
        'index': range(len(plot_data['points'])),
        'detected': plot_data.get('detected_status', [0] * len(plot_data['points']))
    })
    
    # Identify neighbor indices and GUAPOS indices
    neighbor_indices = set()
    guapos_indices = set()
    detected_guapos_indices = set()
    not_detected_guapos_indices = set()
    
    # Get GUAPOS molecules
    guapos_mask = df_points['database'] == 'GUAPOS'
    guapos_indices = set(df_points[guapos_mask]['index'])
    
    # Separate detected and not detected GUAPOS
    detected_guapos_indices = set(df_points[guapos_mask & (df_points['detected'] == 1)]['index'])
    not_detected_guapos_indices = set(df_points[guapos_mask & (df_points['detected'] == 0)]['index'])
    
    # Get neighbor indices from KNN connections
    if plot_data.get('knn_connections'):
        for connection in plot_data['knn_connections']:
            neighbor_indices.add(connection['neighbor_index'])
    
    # Create the main scatter plot
    fig = go.Figure()
    
    # Add points for each database (excluding GUAPOS and neighbors which will be highlighted separately)
    databases = [db for db in df_points['database'].unique() if db != 'GUAPOS']
    for db in databases:
        db_mask = ((df_points['database'] == db) & 
                  (~df_points['index'].isin(neighbor_indices)) &
                  (~df_points['index'].isin(guapos_indices)))
        db_data = df_points[db_mask]
        
        if len(db_data) > 0:
            # Create hover text
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
                    size=8 if db == 'PubChem' else 12,
                    opacity=0.3 if db == 'PubChem' else 0.7
                ),
                name=db,
                hovertext=hover_text,
                hoverinfo='text',
                showlegend=True
            ))
    
    # Add detected GUAPOS molecules (blue)
    if detected_guapos_indices:
        detected_guapos_data = df_points[df_points['index'].isin(detected_guapos_indices)]
        
        if len(detected_guapos_data) > 0:
            # Create hover text for detected GUAPOS
            guapos_hover_text = []
            for _, row in detected_guapos_data.iterrows():
                text = f"<span style='color: blue; font-weight: bold;'>GUAPOS (DETECTED)</span><br>"
                if row['label']:
                    text += f"Name: {row['label']}<br>"
                if row['formula']:
                    formatted_formula = format_chemical_formula(row['formula'])
                    text += f"Formula: {formatted_formula}"
                guapos_hover_text.append(text)
            
            fig.add_trace(go.Scatter(
                x=detected_guapos_data['x'],
                y=detected_guapos_data['y'],
                mode='markers',
                marker=dict(
                    color='blue',
                    size=15,
                    opacity=0.9,
                    line=dict(color='darkblue', width=2)
                ),
                name='GUAPOS (Detected)',
                hovertext=guapos_hover_text,
                hoverinfo='text',
                showlegend=True
            ))
    
    # Add not detected GUAPOS molecules (gray)
    if not_detected_guapos_indices:
        not_detected_guapos_data = df_points[df_points['index'].isin(not_detected_guapos_indices)]
        
        if len(not_detected_guapos_data) > 0:
            # Create hover text for not detected GUAPOS
            guapos_hover_text = []
            for _, row in not_detected_guapos_data.iterrows():
                text = f"<span style='color: gray; font-weight: bold;'>GUAPOS (NOT DETECTED)</span><br>"
                if row['label']:
                    text += f"Name: {row['label']}<br>"
                if row['formula']:
                    formatted_formula = format_chemical_formula(row['formula'])
                    text += f"Formula: {formatted_formula}"
                guapos_hover_text.append(text)
            
            fig.add_trace(go.Scatter(
                x=not_detected_guapos_data['x'],
                y=not_detected_guapos_data['y'],
                mode='markers',
                marker=dict(
                    color='gray',
                    size=12,
                    opacity=0.7,
                    line=dict(color='darkgray', width=1.5)
                ),
                name='GUAPOS (Not Detected)',
                hovertext=guapos_hover_text,
                hoverinfo='text',
                showlegend=True
            ))
    
    # Highlight neighbor points in red
    if neighbor_indices:
        neighbor_data = df_points[df_points['index'].isin(neighbor_indices)]
        
        if len(neighbor_data) > 0:
            # Create hover text for neighbors
            neighbor_hover_text = []
            for _, row in neighbor_data.iterrows():
                text = f"<span style='color: red; font-weight: bold;'>NEIGHBOR CANDIDATE</span><br>"
                text += f"Database: {row['database']}<br>"
                if row['label']:
                    text += f"Name: {row['label']}<br>"
                if row['formula']:
                    formatted_formula = format_chemical_formula(row['formula'])
                    text += f"Formula: {formatted_formula}"
                neighbor_hover_text.append(text)
            
            fig.add_trace(go.Scatter(
                x=neighbor_data['x'],
                y=neighbor_data['y'],
                mode='markers',
                marker=dict(
                    color='red',
                    size=15,
                    opacity=0.9,
                    line=dict(color='darkred', width=2)
                ),
                name='Neighbor Candidates',
                hovertext=neighbor_hover_text,
                hoverinfo='text',
                showlegend=True
            ))
    
    # Add KNN connections if they exist
    if plot_data.get('knn_connections'):
        connections = plot_data['knn_connections']
        for connection in connections:
            guapos_idx = connection['guapos_index']
            neighbor_idx = connection['neighbor_index']
            
            if (guapos_idx < len(plot_data['points']) and 
                neighbor_idx < len(plot_data['points'])):
                
                x0, y0 = plot_data['points'][guapos_idx]
                x1, y1 = plot_data['points'][neighbor_idx]
                
                fig.add_trace(go.Scatter(
                    x=[x0, x1],
                    y=[y0, y1],
                    mode='lines',
                    line=dict(color='gray', width=1.5, dash='dash'),
                    opacity=0.7,
                    showlegend=False,
                    hoverinfo='skip'
                ))
    
    # Update layout
    fig.update_layout(
        title="Chemical Space Visualization - PCA/UMAP Projection",
        xaxis_title="Dimension 1",
        yaxis_title="Dimension 2",
        hovermode='closest',
        showlegend=True,
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        ),
        height=700
    )
    
    fig.update_xaxes(showgrid=False, zeroline=False, showticklabels=False)
    fig.update_yaxes(showgrid=False, zeroline=False, showticklabels=False)
    
    return fig

def create_database_summary(plot_data):
    """Create database summary statistics"""
    databases = plot_data['databases']
    db_counts = pd.Series(databases).value_counts().reset_index()
    db_counts.columns = ['Database', 'Count']
    
    fig = px.bar(
        db_counts,
        x='Database',
        y='Count',
        color='Database',
        title="Molecule Distribution by Database",
        color_discrete_map={
            'GUAPOS': '#1f77b4',
            'TMC_1': '#ff7f0e',
            'All_Discoveries': '#2ca02c',
            'KIDA': '#FFFF00',
            'PubChem': '#9467bd'
        }
    )
    
    fig.update_layout(
        xaxis_title="Database",
        yaxis_title="Number of Molecules",
        showlegend=False
    )
    
    return fig

def main():
    st.markdown('<h1 class="main-header">🧪 Chemical Space Visualization</h1>', unsafe_allow_html=True)
    
    # Sidebar
    with st.sidebar:
        st.header("📊 Data Configuration")
        
        # File upload
        uploaded_file = st.file_uploader(
            "Upload plot data JSON",
            type=['json'],
            help="Upload the plot_data.json file generated by the analysis script"
        )
        
        use_sample_data = st.checkbox("Use sample data for demonstration", value=False)
        
        if use_sample_data:
            plot_data = create_sample_data()
            st.success("✅ Using sample data for demonstration!")
        elif uploaded_file is not None:
            try:
                # Save uploaded file temporarily
                with open("temp_plot_data.json", "wb") as f:
                    f.write(uploaded_file.getvalue())
                
                # Load with validation
                plot_data, message = load_plot_data("temp_plot_data.json")
                if plot_data:
                    st.success(f"✅ {message}")
                else:
                    st.error(f"❌ {message}")
                    plot_data = None
            except Exception as e:
                st.error(f"❌ Error loading file: {str(e)}")
                plot_data = None
        else:
            # Try to load from default path
            try:
                plot_data, message = load_plot_data("plot_data.json")
                if plot_data:
                    st.info(f"📁 {message}")
                else:
                    st.warning(f"⚠️ {message}")
            except Exception as e:
                st.error(f"❌ Error: {str(e)}")
                plot_data = None
        
        if plot_data:
            st.markdown("---")
            st.header("📋 Data Summary")
            
            total_molecules = len(plot_data['points'])
            st.write(f"**Total Molecules:** {total_molecules:,}")
            
            # Database distribution
            db_counts = pd.Series(plot_data['databases']).value_counts()
            st.write("**Database Distribution:**")
            for db, count in db_counts.items():
                st.write(f"• {db}: {count} molecules")
            
            # GUAPOS detection status
            if 'detected_status' in plot_data:
                guapos_indices = [i for i, db in enumerate(plot_data['databases']) if db == 'GUAPOS']
                if guapos_indices:
                    detected_count = sum(plot_data['detected_status'][i] for i in guapos_indices)
                    st.write(f"• GUAPOS Detected: {detected_count}/{len(guapos_indices)}")
            
            # KNN information
            if plot_data.get('knn_connections'):
                st.write(f"• KNN Connections: {len(plot_data['knn_connections'])}")
            
            st.markdown("---")
            st.header("🔍 Interaction Guide")
            st.info("""
            - **Hover** over points to see molecule details
            - **Blue points** = Detected GUAPOS molecules
            - **Gray points** = Not detected GUAPOS molecules
            - **Red points** = Neighbor candidates
            - **Click** on points to view detailed information
            - **Zoom** with mouse wheel or touchpad
            - **Pan** by dragging the plot
            """)
    
    # Main content
    if plot_data is None:
        st.warning("""
        ## Welcome to Chemical Space Visualization!
        
        To get started:
        
        1. **Upload** a plot_data.json file using the sidebar
        2. **Or** check "Use sample data" to see a demonstration
        3. **Or** place plot_data.json in the same directory
        
        ### Common Issues:
        - JSON file might be corrupt or too large
        - Try using the sample data first to verify the app works
        - Check that your JSON file has the correct format
        
        ### Required JSON structure:
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
        
        # Debug information
        st.markdown("---")
        st.header("🛠️ Debug Information")
        
        if uploaded_file:
            st.write("**Uploaded file info:**")
            st.write(f"Name: {uploaded_file.name}")
            st.write(f"Size: {uploaded_file.size} bytes")
            st.write(f"Type: {uploaded_file.type}")
        
        # Show sample data structure
        st.write("**Sample data structure:**")
        st.json({
            "points": [[0.1, 0.2], [0.3, 0.4]],
            "databases": ["GUAPOS", "PubChem"],
            "labels": ["Molecule1", "Molecule2"],
            "formulas": ["H2O", "C6H12O6"],
            "colors": ["#1f77b4", "#9467bd"]
        })
        
        return
    
    # Create tabs for different views
    tab1, tab2 = st.tabs(["📈 Interactive Map", "📊 Statistics"])
    
    with tab1:
        st.markdown('<h2 class="sub-header">Interactive Chemical Space Map</h2>', unsafe_allow_html=True)
        
        # Create interactive plot
        fig = create_interactive_plot(plot_data)
        
        # Display the plot
        st.plotly_chart(fig, use_container_width=True)
        
        # Add information about the visualization
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("""
            <div class="info-box">
            <h4>About this Visualization</h4>
            <p>This map shows molecules projected into 2D space using dimensionality reduction. 
            Each point represents a molecule, colored by its source database.</p>
            </div>
            """, unsafe_allow_html=True)
        
        with col2:
            st.markdown("""
            <div class="info-box">
            <h4>Legend</h4>
            <p>• <strong style='color: blue;'>Blue</strong>: Detected GUAPOS molecules<br>
            • <strong style='color: gray;'>Gray</strong>: Not detected GUAPOS<br>
            • <strong style='color: red;'>Red</strong>: Neighbor candidates<br>
            • <strong>Other colors</strong>: Different databases<br>
            • <strong>Dashed lines</strong>: KNN connections</p>
            </div>
            """, unsafe_allow_html=True)
    
    with tab2:
        st.markdown('<h2 class="sub-header">Database Statistics</h2>', unsafe_allow_html=True)
        
        # Create statistics visualizations
        col1, col2 = st.columns(2)
        
        with col1:
            # Database distribution
            db_fig = create_database_summary(plot_data)
            st.plotly_chart(db_fig, use_container_width=True)
        
        with col2:
            # Statistics cards
            total_molecules = len(plot_data['points'])
            unique_dbs = len(set(plot_data['databases']))
            
            st.metric("Total Molecules", f"{total_molecules:,}")
            st.metric("Unique Databases", unique_dbs)
            
            # GUAPOS statistics
            guapos_count = plot_data['databases'].count('GUAPOS')
            if guapos_count > 0:
                detected = 0
                if 'detected_status' in plot_data:
                    guapos_indices = [i for i, db in enumerate(plot_data['databases']) if db == 'GUAPOS']
                    detected = sum(plot_data['detected_status'][i] for i in guapos_indices)
                st.metric("GUAPOS Molecules", guapos_count)
                st.metric("Detected GUAPOS", f"{detected}/{guapos_count}")
            
            # KNN statistics
            if plot_data.get('knn_connections'):
                st.metric("KNN Connections", len(plot_data['knn_connections']))

if __name__ == "__main__":
    main()
