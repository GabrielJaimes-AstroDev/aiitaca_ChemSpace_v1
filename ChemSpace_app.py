import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import json
import pickle
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
</style>
""", unsafe_allow_html=True)

def find_plot_files():
    """Automatically find plot data files in the current directory"""
    plot_files = []
    for file in os.listdir('.'):
        if file.endswith('.json') or file.endswith('.pkl'):
            plot_files.append(file)
    return sorted(plot_files)

def load_plot_data(file_path):
    """Load plot data from JSON or pickle file"""
    try:
        if file_path.endswith('.json'):
            with open(file_path, 'r') as f:
                return json.load(f)
        elif file_path.endswith('.pkl'):
            with open(file_path, 'rb') as f:
                return pickle.load(f)
    except Exception as e:
        st.error(f"Error loading file {file_path}: {str(e)}")
    return None

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
    """Create interactive Plotly visualization with consistent point sizes"""
    
    fig = go.Figure()
    
    # Group points by database for better organization
    points_by_database = {}
    for point in plot_data.get('points', []):
        db = point.get('database', 'Unknown')
        if db not in points_by_database:
            points_by_database[db] = []
        points_by_database[db].append(point)
    
    # Color mapping
    color_map = plot_data.get('color_scheme', {
        'GUAPOS': '#1f77b4',
        'TMC_1': '#ff7f0e',
        'All_Discoveries': '#2ca02c',
        'KIDA': '#FFFF00',
        'PubChem': '#9467bd'
    })
    
    # Add traces for each database with consistent point sizes
    for db, points in points_by_database.items():
        if not points:
            continue
            
        x = [p['x'] for p in points]
        y = [p['y'] for p in points]
        
        # Set consistent point sizes (smaller than before)
        base_size = 6
        if db == 'GUAPOS':
            # Smaller size for GUAPOS points
            sizes = [base_size * 1.2 for _ in points]  # Slightly larger but not huge
        elif db == plot_data.get('neighbor_database', 'PubChem'):
            # Smaller size for neighbor points
            sizes = [base_size for _ in points]
        else:
            sizes = [base_size for _ in points]
        
        # Create hover text
        hover_text = []
        for point in points:
            text = f"<b>{point.get('name', 'Unknown')}</b><br>"
            if point.get('formula'):
                text += f"Formula: {point['formula']}<br>"
            text += f"Database: {point.get('database', 'Unknown')}<br>"
            if point.get('detected') is not None:
                status = "Detected" if point['detected'] else "Not Detected"
                text += f"Status: {status}<br>"
            if point.get('is_neighbor'):
                text += "Type: KNN Neighbor<br>"
            hover_text.append(text)
        
        # Get color from mapping or use default
        color = color_map.get(db, '#cccccc')
        
        fig.add_trace(go.Scatter(
            x=x,
            y=y,
            mode='markers',
            marker=dict(
                size=sizes,
                color=color,
                opacity=0.7,
                line=dict(width=0.5, color='DarkSlateGrey')
            ),
            name=db,
            text=hover_text,
            hoverinfo='text',
            hovertemplate="%{text}<extra></extra>"
        ))
    
    # Add KNN connections if they exist
    if 'connections' in plot_data:
        connections = plot_data['connections']
        for connection in connections:
            # Find the actual points
            guapos_point = None
            neighbor_point = None
            
            for point in plot_data['points']:
                if (point.get('database') == 'GUAPOS' and 
                    point.get('detected') and 
                    'index' in point and point['index'] == connection.get('guapos_index')):
                    guapos_point = point
                elif ('index' in point and point['index'] == connection.get('neighbor_index')):
                    neighbor_point = point
            
            if guapos_point and neighbor_point:
                fig.add_trace(go.Scatter(
                    x=[guapos_point['x'], neighbor_point['x']],
                    y=[guapos_point['y'], neighbor_point['y']],
                    mode='lines',
                    line=dict(color='gray', width=1, dash='dash'),
                    opacity=0.6,
                    showlegend=False,
                    hoverinfo='none'
                ))
    
    # Add labels if they exist
    if 'labels' in plot_data:
        for label in plot_data['labels']:
            fig.add_annotation(
                x=label['x'],
                y=label['y'],
                text=label['text'],
                showarrow=False,
                font=dict(size=label.get('fontsize', 8), color='black'),
                opacity=label.get('alpha', 0.8),
                bgcolor='rgba(255,255,255,0.7)',
                bordercolor='black',
                borderwidth=1,
                borderpad=2
            )
    
    # Update layout
    title = plot_data.get('title', 'Chemical Space Analysis')
    fig.update_layout(
        title=title,
        showlegend=True,
        legend=dict(title='Databases', itemsizing='constant'),
        hovermode='closest',
        width=1000,
        height=700,
        plot_bgcolor='white'
    )
    
    # Remove axis labels and ticks
    fig.update_xaxes(showticklabels=False, showgrid=False, zeroline=False)
    fig.update_yaxes(showticklabels=False, showgrid=False, zeroline=False)
    
    return fig

def create_database_summary(plot_data):
    """Create database summary statistics"""
    db_counts = {}
    for point in plot_data.get('points', []):
        db = point.get('database', 'Unknown')
        db_counts[db] = db_counts.get(db, 0) + 1
    
    summary_df = pd.DataFrame({
        'Database': list(db_counts.keys()),
        'Count': list(db_counts.values())
    }).sort_values('Count', ascending=False)
    
    # Create bar chart with consistent colors
    color_map = plot_data.get('color_scheme', {
        'GUAPOS': '#1f77b4',
        'TMC_1': '#ff7f0e',
        'All_Discoveries': '#2ca02c',
        'KIDA': '#FFFF00',
        'PubChem': '#9467bd'
    })
    
    fig = px.bar(
        summary_df,
        x='Database',
        y='Count',
        color='Database',
        title="Molecule Distribution by Database",
        color_discrete_map=color_map
    )
    
    fig.update_layout(
        xaxis_title="Database",
        yaxis_title="Number of Molecules",
        showlegend=False
    )
    
    return fig, summary_df

def create_knn_results_table(plot_data):
    """Create a table of KNN results"""
    if 'knn_results' not in plot_data:
        return None
    
    results = []
    for result in plot_data['knn_results']:
        if result.get('detected', 0) == 1:
            for neighbor in result.get('neighbors', []):
                results.append({
                    'GUAPOS Molecule': result['guapos_name'],
                    'GUAPOS Formula': result['guapos_formula'],
                    'Neighbor Name': neighbor['name'],
                    'Neighbor Formula': neighbor['formula'],
                    'Neighbor Database': neighbor['database'],
                    'Distance': f"{neighbor['distance']:.4f}",
                    'Rank': neighbor['rank']
                })
    
    if results:
        return pd.DataFrame(results)
    return None

def main():
    st.markdown('<h1 class="main-header">🧪 Chemical Space Visualization</h1>', unsafe_allow_html=True)
    
    # Automatically find plot files
    plot_files = find_plot_files()
    
    # Sidebar
    with st.sidebar:
        st.header("📊 Data Selection")
        
        if not plot_files:
            st.warning("No plot data files found. Please add JSON or PKL files to the repository.")
            return
        
        # File selection
        selected_file = st.selectbox(
            "Select plot data file:",
            options=plot_files,
            help="Select the plot data file to visualize"
        )
        
        if selected_file:
            plot_data = load_plot_data(selected_file)
            
            if plot_data:
                st.success(f"✅ Loaded: {selected_file}")
                
                st.markdown("---")
                st.header("🎨 Display Options")
                
                # Database filter
                all_dbs = list(set(point.get('database', 'Unknown') 
                                 for point in plot_data.get('points', [])))
                selected_dbs = st.multiselect(
                    "Filter databases:",
                    options=all_dbs,
                    default=all_dbs,
                    help="Select which databases to display"
                )
                
                # Point size adjustment
                point_size = st.slider(
                    "Point size:",
                    min_value=4,
                    max_value=12,
                    value=6,
                    step=1,
                    help="Adjust the size of all points"
                )
                
                # Show/hide connections
                show_connections = st.checkbox(
                    "Show KNN connections",
                    value=True,
                    help="Display lines between GUAPOS molecules and their neighbors"
                )
                
                st.markdown("---")
                st.header("📊 Dataset Info")
                
                # Show basic stats
                total_points = len(plot_data.get('points', []))
                st.write(f"**Total molecules:** {total_points}")
                
                if 'knn_results' in plot_data:
                    detected = sum(1 for r in plot_data['knn_results'] if r.get('detected', 0) == 1)
                    st.write(f"**Detected GUAPOS:** {detected}")
                
                neighbor_db = plot_data.get('neighbor_database', 'Unknown')
                st.write(f"**Neighbor database:** {neighbor_db}")
    
    # Main content
    if 'plot_data' not in locals():
        st.info("👈 Please select a plot data file from the sidebar")
        
        if plot_files:
            st.write("**Available files:**")
            for file in plot_files:
                st.write(f"• `{file}`")
        
        st.markdown("""
        ### 📋 How to Use
        
        1. **Add plot data files** (JSON or PKL) to your GitHub repository
        2. **Select a file** from the sidebar dropdown
        3. **Explore** the interactive visualization
        4. **Adjust settings** using the sidebar controls
        
        ### 🧪 Expected Data Format
        
        The plot data files should contain:
        - Molecular coordinates and properties
        - Database information
        - KNN analysis results
        - Visualization parameters
        """)
        return
    
    # Create tabs for different views
    tab1, tab2, tab3, tab4 = st.tabs(["📈 Interactive Map", "📊 Statistics", "🔍 KNN Results", "ℹ️ About"])
    
    with tab1:
        st.markdown('<h2 class="sub-header">Interactive Chemical Space Map</h2>', unsafe_allow_html=True)
        
        # Create interactive plot
        fig = create_interactive_plot(plot_data)
        
        # Display the plot
        st.plotly_chart(fig, use_container_width=True)
        
        # Add context information
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("""
            <div class="info-box">
            <h4>About this Visualization</h4>
            <p>This map shows molecules projected into 2D space using PCA and UMAP. 
            Colors represent different databases, and connections show molecular similarities.</p>
            </div>
            """, unsafe_allow_html=True)
        
        with col2:
            st.markdown("""
            <div class="info-box">
            <h4>How to Interpret</h4>
            <p>• <strong>Closer points</strong> = More similar molecules<br>
            • <strong>Colors</strong> = Source database<br>
            • <strong>Dashed lines</strong> = KNN similarity connections<br>
            • <strong>Hover</strong> for molecule details</p>
            </div>
            """, unsafe_allow_html=True)
    
    with tab2:
        st.markdown('<h2 class="sub-header">Database Statistics</h2>', unsafe_allow_html=True)
        
        # Create statistics visualizations
        db_fig, summary_df = create_database_summary(plot_data)
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.plotly_chart(db_fig, use_container_width=True)
        
        with col2:
            st.dataframe(
                summary_df,
                use_container_width=True,
                hide_index=True
            )
            
            # Additional stats
            total_molecules = len(plot_data.get('points', []))
            unique_dbs = len(summary_df)
            
            st.metric("Total Molecules", total_molecules)
            st.metric("Unique Databases", unique_dbs)
    
    with tab3:
        st.markdown('<h2 class="sub-header">KNN Analysis Results</h2>', unsafe_allow_html=True)
        
        knn_table = create_knn_results_table(plot_data)
        if knn_table is not None:
            st.dataframe(
                knn_table,
                use_container_width=True,
                hide_index=True,
                height=400
            )
            
            # Download button
            csv = knn_table.to_csv(index=False)
            st.download_button(
                "📥 Download KNN Results",
                data=csv,
                file_name="knn_results.csv",
                mime="text/csv"
            )
        else:
            st.info("No KNN results found in the selected file.")
    
    with tab4:
        st.markdown('<h2 class="sub-header">About This Visualization</h2>', unsafe_allow_html=True)
        
        st.markdown("""
        ### 🧪 Chemical Space Analysis
        
        This interactive visualization shows the results of chemical space analysis using:
        
        - **PCA (Principal Component Analysis)**: Linear dimensionality reduction
        - **UMAP (Uniform Manifold Approximation and Projection)**: Non-linear dimensionality reduction
        - **KNN (K-Nearest Neighbors)**: Similarity search algorithm
        
        ### 🎨 Color Coding
        
        - **Blue**: GUAPOS molecules
        - **Orange**: TMC_1 database
        - **Green**: All_Discoveries database
        - **Yellow**: KIDA database
        - **Purple**: PubChem database
        
        ### 🔍 Interactive Features
        
        - **Hover** over points to see molecule details
        - **Zoom** and **pan** to explore different regions
        - **Filter** databases using the sidebar
        - **Download** results for further analysis
        
        ### 📁 Data Sources
        
        The visualization automatically detects plot data files in the repository.
        Supported formats: JSON and PKL files generated by the analysis pipeline.
        """)

if __name__ == "__main__":
    main()
