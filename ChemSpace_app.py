import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import json
import pickle
import re

# Page configuration
st.set_page_config(
    page_title="Molecular Space Analysis",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS
st.markdown("""
<style>
    .main-header {
        font-size: 3rem;
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
        background-color: #f0f2f6;
        padding: 1rem;
        border-radius: 0.5rem;
        margin-bottom: 1rem;
    }
    .molecule-card {
        border: 1px solid #ddd;
        border-radius: 0.5rem;
        padding: 1rem;
        margin-bottom: 1rem;
        background-color: white;
    }
</style>
""", unsafe_allow_html=True)

def load_data(file_path):
    """Load data from JSON or pickle file"""
    if file_path.endswith('.json'):
        with open(file_path, 'r') as f:
            return json.load(f)
    elif file_path.endswith('.pkl'):
        with open(file_path, 'rb') as f:
            return pickle.load(f)
    return None

def reformat_chemical_formula(formula):
    """Format chemical formula with subscripts"""
    try:
        if not formula or pd.isna(formula) or formula.lower() == 'nan':
            return ""
        formula = str(formula).strip()
        formula = re.sub(r'(\d+)', r'<sub>\1</sub>', formula)
        formula = re.sub(r'([+-]\d*)', r'<sup>\1</sup>', formula)
        return formula
    except:
        return formula

def create_interactive_plot(data):
    """Create interactive Plotly visualization"""
    molecules_df = pd.DataFrame(data['molecules'])
    
    # Create figure
    fig = go.Figure()
    
    # Group by database and add traces
    databases = molecules_df['database'].unique()
    colors = data.get('visualization_settings', {}).get('db_colors', {})
    
    for db in databases:
        db_data = molecules_df[molecules_df['database'] == db]
        color = colors.get(db, '#cccccc')
        
        fig.add_trace(go.Scatter(
            x=db_data['x'],
            y=db_data['y'],
            mode='markers',
            name=db,
            marker=dict(
                size=8 if db != 'PubChem' else 6,
                color=color,
                opacity=0.7 if db != 'PubChem' else 0.3,
                line=dict(width=0.5, color='black')
            ),
            text=db_data.apply(lambda row: f"""
                <b>Name:</b> {row['name']}<br>
                <b>Formula:</b> {reformat_chemical_formula(row['formula'])}<br>
                <b>Database:</b> {row['database']}<br>
                <b>Detected:</b> {row.get('detected', 'N/A')}
            """, axis=1),
            hoverinfo='text',
            customdata=db_data[['name', 'formula', 'database', 'detected']]
        ))
    
    # Add GUAPOS molecules with special styling
    if 'GUAPOS' in databases:
        guapos_data = molecules_df[molecules_df['database'] == 'GUAPOS']
        detected_data = guapos_data[guapos_data['detected'] == 1]
        not_detected_data = guapos_data[guapos_data['detected'] == 0]
        
        if not detected_data.empty:
            fig.add_trace(go.Scatter(
                x=detected_data['x'],
                y=detected_data['y'],
                mode='markers+text',
                name='GUAPOS (Detected)',
                marker=dict(
                    size=12,
                    color='blue',
                    opacity=0.9,
                    line=dict(width=2, color='black')
                ),
                text=[reformat_chemical_formula(f) for f in detected_data['formula']],
                textposition='top center',
                textfont=dict(size=10, color='black'),
                hoverinfo='text',
                customdata=detected_data[['name', 'formula', 'database', 'detected']]
            ))
        
        if not not_detected_data.empty:
            fig.add_trace(go.Scatter(
                x=not_detected_data['x'],
                y=not_detected_data['y'],
                mode='markers',
                name='GUAPOS (Not Detected)',
                marker=dict(
                    size=10,
                    color='gray',
                    opacity=0.6,
                    line=dict(width=1, color='black')
                ),
                hoverinfo='text',
                customdata=not_detected_data[['name', 'formula', 'database', 'detected']]
            ))
    
    # Update layout
    fig.update_layout(
        title=f"Molecular Space Analysis - {data['metadata']['neighbor_database']} Neighbors",
        xaxis_title="UMAP Dimension 1",
        yaxis_title="UMAP Dimension 2",
        showlegend=True,
        hovermode='closest',
        width=1000,
        height=700,
        plot_bgcolor='white'
    )
    
    fig.update_xaxes(showgrid=False, zeroline=False)
    fig.update_yaxes(showgrid=False, zeroline=False)
    
    return fig

def create_statistics_dashboard(data):
    """Create statistics dashboard"""
    molecules_df = pd.DataFrame(data['molecules'])
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("Total Molecules", len(molecules_df))
    
    with col2:
        guapos_count = len(molecules_df[molecules_df['database'] == 'GUAPOS'])
        st.metric("GUAPOS Molecules", guapos_count)
    
    with col3:
        neighbor_db = data['metadata']['neighbor_database']
        neighbor_count = len(molecules_df[molecules_df['database'] == neighbor_db])
        st.metric(f"{neighbor_db} Molecules", neighbor_count)
    
    # Database distribution
    st.subheader("Database Distribution")
    db_counts = molecules_df['database'].value_counts()
    fig_db = px.pie(
        values=db_counts.values,
        names=db_counts.index,
        title="Molecules by Database"
    )
    st.plotly_chart(fig_db, use_container_width=True)
    
    # Detection status for GUAPOS
    if 'GUAPOS' in molecules_df['database'].unique():
        st.subheader("GUAPOS Detection Status")
        guapos_data = molecules_df[molecules_df['database'] == 'GUAPOS']
        detected_counts = guapos_data['detected'].value_counts()
        
        col1, col2 = st.columns(2)
        
        with col1:
            fig_detection = px.pie(
                values=detected_counts.values,
                names=detected_counts.index.map({1: 'Detected', 0: 'Not Detected'}),
                title="GUAPOS Detection Status"
            )
            st.plotly_chart(fig_detection, use_container_width=True)
        
        with col2:
            if data.get('knn_results'):
                detected_with_neighbors = sum(1 for res in data['knn_results'] 
                                           if res.get('detected') == 1 and res.get('neighbors'))
                st.metric("GUAPOS with Neighbors", detected_with_neighbors)

def display_neighbors_table(data):
    """Display neighbors table"""
    if not data.get('knn_results'):
        st.warning("No KNN results available")
        return
    
    knn_results = data['knn_results']
    neighbor_db = data['metadata']['neighbor_database']
    
    st.subheader(f"KNN Analysis Results - {neighbor_db} Neighbors")
    
    # Filter only detected molecules with neighbors
    detected_results = [res for res in knn_results if res.get('detected') == 1 and res.get('neighbors')]
    
    if not detected_results:
        st.info("No detected molecules with neighbors found")
        return
    
    # Create expandable sections for each GUAPOS molecule
    for i, result in enumerate(detected_results):
        with st.expander(f"{result['guapos_name']} - {reformat_chemical_formula(result['guapos_formula'])}"):
            col1, col2 = st.columns(2)
            
            with col1:
                st.markdown(f"**SMILES:** `{result['guapos_smiles']}`")
                st.markdown(f"**Status:** Detected")
            
            with col2:
                if result['neighbors']:
                    st.markdown(f"**Nearest Neighbor:** {result['neighbors'][0]['name']}")
                    st.markdown(f"**Distance:** {result['neighbors'][0]['distance']:.4f}")
            
            # Display neighbors table
            neighbors_df = pd.DataFrame(result['neighbors'])
            st.dataframe(
                neighbors_df[['rank', 'name', 'formula', 'distance', 'database']],
                use_container_width=True
            )

def main():
    # Header
    st.markdown('<h1 class="main-header">🧪 Molecular Space Analysis</h1>', unsafe_allow_html=True)
    
    # Sidebar for file upload
    with st.sidebar:
        st.header("Data Input")
        uploaded_file = st.file_uploader(
            "Upload analysis data",
            type=['json', 'pkl'],
            help="Upload the JSON or pickle file generated by the analysis"
        )
        
        if uploaded_file:
            # Save uploaded file temporarily
            with open("temp_data.file", "wb") as f:
                f.write(uploaded_file.getbuffer())
            
            # Load data
            try:
                data = load_data("temp_data.file")
                st.success("Data loaded successfully!")
                
                # Display metadata
                st.subheader("Dataset Information")
                st.write(f"**Total molecules:** {data['metadata']['total_molecules']}")
                st.write(f"**Neighbor database:** {data['metadata']['neighbor_database']}")
                st.write(f"**Analysis date:** {data['metadata']['timestamp']}")
                
                # Database counts
                st.subheader("Database Counts")
                for db, count in data['metadata']['databases_present'].items():
                    st.write(f"{db}: {count}")
                    
            except Exception as e:
                st.error(f"Error loading data: {str(e)}")
                data = None
        else:
            data = None
            st.info("Please upload a analysis data file to begin")
    
    # Main content
    if data:
        # Create tabs
        tab1, tab2, tab3 = st.tabs(["Interactive Plot", "Statistics", "Neighbors Table"])
        
        with tab1:
            st.subheader("Interactive Molecular Space Visualization")
            fig = create_interactive_plot(data)
            st.plotly_chart(fig, use_container_width=True)
            
            # Additional controls
            col1, col2 = st.columns(2)
            
            with col1:
                st.subheader("Plot Controls")
                show_labels = st.checkbox("Show molecule labels", value=True)
                opacity = st.slider("Marker opacity", 0.1, 1.0, 0.7)
            
            with col2:
                st.subheader("Database Filter")
                databases = list(set([m['database'] for m in data['molecules']]))
                selected_dbs = st.multiselect(
                    "Select databases to show",
                    databases,
                    default=databases
                )
        
        with tab2:
            create_statistics_dashboard(data)
        
        with tab3:
            display_neighbors_table(data)
            
    else:
        # Welcome screen with instructions
        st.markdown("""
        <div class="info-box">
            <h3>Welcome to Molecular Space Analysis</h3>
            <p>This application visualizes molecular space analysis data including:</p>
            <ul>
                <li>Interactive 2D molecular space visualization</li>
                <li>Statistical analysis of molecular distributions</li>
                <li>KNN neighbor results tables</li>
                <li>Database-specific molecular information</li>
            </ul>
            <p><b>To get started:</b></p>
            <ol>
                <li>Run the molecular analysis script to generate data</li>
                <li>Upload the resulting JSON or pickle file using the sidebar</li>
                <li>Explore the results through the interactive tabs</li>
            </ol>
        </div>
        """, unsafe_allow_html=True)
        
        # Example of what the interface will show
        col1, col2 = st.columns(2)
        
        with col1:
            st.subheader("Features")
            st.markdown("""
            - **Interactive Plot**: Zoom, pan, and hover over molecules
            - **Color-coded databases**: Different colors for each molecular database
            - **Detection status**: Visual distinction between detected and non-detected molecules
            - **Detailed tooltips**: Comprehensive information on hover
            """)
        
        with col2:
            st.subheader("Supported Data")
            st.markdown("""
            - Multiple molecular databases (GUAPOS, PubChem, TMC, etc.)
            - Detection status information
            - KNN analysis results
            - Molecular fingerprints and coordinates
            - Chemical formulas and SMILES strings
            """)

if __name__ == "__main__":
    main()
