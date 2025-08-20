import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import json
import pickle
import re
from io import StringIO

def load_plot_data(file_path):
    """Load plot data from JSON or pickle file"""
    if file_path.endswith('.json'):
        with open(file_path, 'r') as f:
            return json.load(f)
    elif file_path.endswith('.pkl'):
        with open(file_path, 'rb') as f:
            return pickle.load(f)
    return None

def parse_chemical_formula(formula_text):
    """Parse chemical formula with subscripts and superscripts"""
    if not formula_text:
        return ""
    
    # Convert LaTeX-style formatting to HTML
    formula = formula_text.replace('$_{', '<sub>').replace('}$', '</sub>')
    formula = formula.replace('$^{', '<sup>').replace('}$', '</sup>')
    formula = formula.replace('$', '')
    
    return formula

def create_interactive_plot(plot_data):
    """Create an interactive Plotly visualization from the plot data"""
    
    # Create figure
    fig = go.Figure()
    
    # Group points by database for better organization
    points_by_database = {}
    for point in plot_data.get('points', []):
        db = point.get('database', 'Unknown')
        if db not in points_by_database:
            points_by_database[db] = []
        points_by_database[db].append(point)
    
    # Add traces for each database
    color_map = plot_data.get('color_scheme', {})
    
    for db, points in points_by_database.items():
        if not points:
            continue
            
        x = [p['x'] for p in points]
        y = [p['y'] for p in points]
        colors = [p.get('color', color_map.get(db, '#cccccc')) for p in points]
        sizes = [p.get('size', 20) for p in points]
        alphas = [p.get('alpha', 0.7) for p in points]
        names = [p.get('name', 'Unknown') for p in points]
        formulas = [p.get('formula', '') for p in points]
        smiles = [p.get('smiles', '') for p in points]
        
        # Create hover text
        hover_text = []
        for i, point in enumerate(points):
            text = f"<b>{point.get('name', 'Unknown')}</b><br>"
            if point.get('formula'):
                text += f"Formula: {point['formula']}<br>"
            if point.get('database'):
                text += f"Database: {point['database']}<br>"
            if point.get('detected') is not None:
                status = "Detected" if point['detected'] else "Not Detected"
                text += f"Status: {status}<br>"
            if point.get('is_neighbor'):
                text += "Type: KNN Neighbor<br>"
            hover_text.append(text)
        
        # Add scatter trace
        fig.add_trace(go.Scatter(
            x=x,
            y=y,
            mode='markers',
            marker=dict(
                size=sizes,
                color=colors,
                opacity=alphas,
                line=dict(width=0.5, color='DarkSlateGrey')
            ),
            name=db,
            text=hover_text,
            hoverinfo='text',
            customdata=list(zip(names, formulas, smiles)),
            hovertemplate=(
                "%{text}<br>" +
                "<extra></extra>"
            )
        ))
    
    # Add KNN connections if they exist
    if 'connections' in plot_data:
        connections = plot_data['connections']
        for connection in connections:
            guapos_idx = connection['guapos_index']
            neighbor_idx = connection['neighbor_index']
            
            # Find the points
            guapos_point = None
            neighbor_point = None
            
            for point in plot_data['points']:
                if point.get('database') == 'GUAPOS' and point.get('detected') and point.get('name') == plot_data['points'][guapos_idx].get('name'):
                    guapos_point = point
                elif point.get('index') == neighbor_idx:
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
    
    # Add labels
    if 'labels' in plot_data:
        for label in plot_data['labels']:
            fig.add_annotation(
                x=label['x'],
                y=label['y'],
                text=label['text'],
                showarrow=False,
                font=dict(
                    size=label.get('fontsize', 8),
                    color='black'
                ),
                opacity=label.get('alpha', 0.8),
                bgcolor='white' if label.get('type') == 'guapos' else 'rgba(255,255,255,0.7)',
                bordercolor='black',
                borderwidth=1 if label.get('type') == 'guapos' else 0,
                borderpad=4 if label.get('type') == 'guapos' else 2
            )
    
    # Update layout
    fig.update_layout(
        title=plot_data.get('title', 'Chemical Space Analysis'),
        showlegend=True,
        legend=dict(
            title='Databases',
            itemsizing='constant'
        ),
        hovermode='closest',
        width=1000,
        height=800,
        plot_bgcolor='white'
    )
    
    # Remove axis labels and ticks
    fig.update_xaxes(showticklabels=False, showgrid=False, zeroline=False)
    fig.update_yaxes(showticklabels=False, showgrid=False, zeroline=False)
    
    return fig

def create_database_summary(plot_data):
    """Create a summary of databases in the plot"""
    db_counts = {}
    for point in plot_data.get('points', []):
        db = point.get('database', 'Unknown')
        db_counts[db] = db_counts.get(db, 0) + 1
    
    summary_df = pd.DataFrame({
        'Database': list(db_counts.keys()),
        'Count': list(db_counts.values())
    }).sort_values('Count', ascending=False)
    
    return summary_df

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
    st.set_page_config(
        page_title="Chemical Space Visualizer",
        page_icon="🧪",
        layout="wide",
        initial_sidebar_state="expanded"
    )
    
    st.title("🧪 Chemical Space Analysis Visualizer")
    st.markdown("""
    Interactive visualization of chemical space analysis using PCA and UMAP projections.
    Upload your exported plot data to explore the molecular landscape.
    """)
    
    # Sidebar for file upload and controls
    with st.sidebar:
        st.header("📁 Data Upload")
        uploaded_file = st.file_uploader(
            "Upload plot data (JSON or PKL)",
            type=['json', 'pkl'],
            help="Upload the exported plot data file from the analysis"
        )
        
        if uploaded_file:
            st.success("File uploaded successfully!")
            
            # Read file content
            file_content = uploaded_file.read()
            
            # Try to load as JSON first, then as pickle
            plot_data = None
            try:
                if uploaded_file.name.endswith('.json'):
                    plot_data = json.loads(file_content.decode('utf-8'))
                else:
                    # For pickle, we need to write to a temporary file
                    import tempfile
                    with tempfile.NamedTemporaryFile(delete=False, suffix='.pkl') as tmp:
                        tmp.write(file_content)
                        tmp.flush()
                        with open(tmp.name, 'rb') as f:
                            plot_data = pickle.load(f)
            except Exception as e:
                st.error(f"Error loading file: {str(e)}")
                plot_data = None
            
            if plot_data:
                st.header("🎛️ Display Options")
                
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
                    "Point size multiplier:",
                    min_value=0.1,
                    max_value=3.0,
                    value=1.0,
                    step=0.1,
                    help="Adjust the size of all points"
                )
                
                # Show/hide connections
                show_connections = st.checkbox(
                    "Show KNN connections",
                    value=True,
                    help="Display lines between GUAPOS molecules and their neighbors"
                )
                
                # Show/hide labels
                show_labels = st.checkbox(
                    "Show molecule labels",
                    value=True,
                    help="Display chemical formulas on the plot"
                )
    
    # Main content area
    if uploaded_file and plot_data:
        # Create tabs for different views
        tab1, tab2, tab3, tab4 = st.tabs(["📊 Visualization", "📈 Statistics", "🔍 KNN Results", "ℹ️ About"])
        
        with tab1:
            st.header("Interactive Chemical Space Map")
            
            # Apply filters if needed
            filtered_plot_data = plot_data.copy()
            if 'points' in filtered_plot_data and selected_dbs:
                filtered_plot_data['points'] = [
                    p for p in filtered_plot_data['points'] 
                    if p.get('database', 'Unknown') in selected_dbs
                ]
            
            # Adjust point sizes
            if 'points' in filtered_plot_data and point_size != 1.0:
                for point in filtered_plot_data['points']:
                    if 'size' in point:
                        point['size'] = point['size'] * point_size
            
            # Create interactive plot
            fig = create_interactive_plot(filtered_plot_data)
            
            # Display the plot
            st.plotly_chart(fig, use_container_width=True)
            
            # Add some controls below the plot
            col1, col2, col3 = st.columns(3)
            with col1:
                st.download_button(
                    "📥 Download Plot as HTML",
                    data=fig.to_html(),
                    file_name="chemical_space_plot.html",
                    mime="text/html"
                )
            with col2:
                st.download_button(
                    "📊 Download Plot as PNG",
                    data=fig.to_image(format="png"),
                    file_name="chemical_space_plot.png",
                    mime="image/png"
                )
        
        with tab2:
            st.header("Dataset Statistics")
            
            # Database summary
            summary_df = create_database_summary(plot_data)
            if not summary_df.empty:
                st.subheader("Database Distribution")
                
                col1, col2 = st.columns(2)
                
                with col1:
                    st.dataframe(
                        summary_df,
                        use_container_width=True,
                        hide_index=True
                    )
                
                with col2:
                    fig_pie = px.pie(
                        summary_df,
                        values='Count',
                        names='Database',
                        title='Database Distribution'
                    )
                    st.plotly_chart(fig_pie, use_container_width=True)
            
            # Additional statistics
            st.subheader("Additional Information")
            
            if 'knn_results' in plot_data:
                detected_count = sum(1 for r in plot_data['knn_results'] if r.get('detected', 0) == 1)
                total_guapos = len(plot_data['knn_results'])
                
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric("Total GUAPOS Molecules", total_guapos)
                with col2:
                    st.metric("Detected Molecules", detected_count)
                with col3:
                    st.metric("Detection Rate", f"{(detected_count/total_guapos*100):.1f}%")
        
        with tab3:
            st.header("KNN Analysis Results")
            
            knn_table = create_knn_results_table(plot_data)
            if knn_table is not None:
                st.dataframe(
                    knn_table,
                    use_container_width=True,
                    hide_index=True
                )
                
                # Download button for KNN results
                csv = knn_table.to_csv(index=False)
                st.download_button(
                    "📥 Download KNN Results as CSV",
                    data=csv,
                    file_name="knn_results.csv",
                    mime="text/csv"
                )
            else:
                st.info("No KNN results found in the uploaded data.")
        
        with tab4:
            st.header("About This Visualization")
            st.markdown("""
            ### 📊 Chemical Space Analysis
            
            This interactive visualization shows the results of chemical space analysis using:
            
            - **PCA (Principal Component Analysis)**: Dimensionality reduction technique
            - **UMAP (Uniform Manifold Approximation and Projection)**: Non-linear dimensionality reduction
            - **KNN (K-Nearest Neighbors)**: Machine learning algorithm for similarity search
            
            ### 🎨 Color Coding
            
            - **Blue**: GUAPOS molecules (Detected)
            - **Gray**: GUAPOS molecules (Not Detected)
            - **Red**: Neighbor molecules from the selected database
            - **Other colors**: Different chemical databases
            
            ### 🔍 Interactive Features
            
            - Hover over points to see molecule details
            - Click and drag to pan around the plot
            - Use the mouse wheel to zoom in/out
            - Toggle databases using the legend
            - Download high-quality images of the plot
            
            ### 📁 Data Export
            
            The visualization uses data exported from the chemical space analysis pipeline,
            including molecular coordinates, properties, and KNN relationships.
            """)
    
    else:
        # Show instructions if no file uploaded
        st.info("👈 Please upload a plot data file to begin visualization")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("""
            ### 📋 How to Use
            
            1. **Run the analysis** using the chemical space analysis code
            2. **Export the plot data** using the export functionality
            3. **Upload the JSON or PKL file** in the sidebar
            4. **Explore** the interactive visualization
            
            ### 📊 Expected File Format
            
            The uploaded file should contain:
            - Molecular coordinates (x, y)
            - Database information
            - Chemical formulas and names
            - KNN relationship data
            - Visualization parameters
            """)
        
        with col2:
            st.markdown("""
            ### 🧪 Supported Databases
            
            - **GUAPOS**: Detected and not detected molecules
            - **PubChem**: Large chemical database
            - **TMC_1**: Theoretical molecular catalog
            - **All_Discoveries**: Experimental discoveries
            - **KIDA**: Kinetic Database for Astrochemistry
            
            ### 🔬 Analysis Features
            
            - Molecular fingerprint generation
            - Dimensionality reduction (PCA + UMAP)
            - Similarity search (KNN)
            - Interactive visualization
            - Data export for external analysis
            """)

if __name__ == "__main__":
    main()
