import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import json
import re
from sklearn.neighbors import NearestNeighbors

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
    .main-title {
        font-size: 1.8rem;
        color: #1f77b4;
        text-align: center;
        margin-bottom: 1rem;
        font-weight: bold;
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

def load_plot_data(file_path):
    """Load plot data from JSON file"""
    try:
        with open(file_path, 'r') as f:
            return json.load(f)
    except Exception as e:
        st.error(f"Error loading plot data: {str(e)}")
        return None

def format_chemical_formula(formula):
    """Format chemical formula with subscripts"""
    if not formula or pd.isna(formula):
        return ""
    
    # Handle subscripts
    formula = re.sub(r'(\d+)', r'<sub>\1</sub>', str(formula))
    
    # Handle ionic charges
    formula = re.sub(r'([+-]\d*)', r'<sup>\1</sup>', formula)
    
    return formula

def find_knn_neighbors(plot_data):
    """Find KNN neighbors from PubChem that are close to GUAPOS points"""
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
    
    # Get GUAPOS points (verificar que existan)
    guapos_mask = df_points['database'] == 'GUAPOS'
    if not guapos_mask.any():
        return []
    
    guapos_points = df_points[guapos_mask][['x', 'y']].values
    guapos_indices = df_points[guapos_mask]['index'].values
    
    # Get PubChem points
    pubchem_mask = (df_points['database'] == 'PubChem')
    if not pubchem_mask.any():
        return []
    
    pubchem_points = df_points[pubchem_mask][['x', 'y']].values
    pubchem_indices = df_points[pubchem_mask]['index'].values
    
    if len(guapos_points) == 0 or len(pubchem_points) == 0:
        return []
    
    # Find nearest neighbors for each GUAPOS point
    n_neighbors = min(10, len(pubchem_points))
    knn = NearestNeighbors(n_neighbors=n_neighbors)
    knn.fit(pubchem_points)
    
    # Get distances and indices of nearest neighbors
    distances, indices = knn.kneighbors(guapos_points)
    
    # Collect unique neighbors (avoid duplicates)
    unique_neighbors = set()
    neighbor_info = []
    
    for i, (guapos_idx, guapos_point) in enumerate(zip(guapos_indices, guapos_points)):
        for j, neighbor_pos in enumerate(indices[i]):
            pubchem_index = pubchem_indices[neighbor_pos]
            
            # Verificar que el índice sea válido
            if (pubchem_index < len(df_points) and 
                pubchem_index not in unique_neighbors):
                
                unique_neighbors.add(pubchem_index)
                
                # Obtener información con validación
                guapos_label = df_points.iloc[guapos_idx]['label'] if guapos_idx < len(df_points) else f"Index {guapos_idx}"
                neighbor_label = df_points.iloc[pubchem_index]['label'] if pubchem_index < len(df_points) else f"Index {pubchem_index}"
                guapos_formula = df_points.iloc[guapos_idx]['formula'] if guapos_idx < len(df_points) else "Unknown"
                neighbor_formula = df_points.iloc[pubchem_index]['formula'] if pubchem_index < len(df_points) else "Unknown"
                
                neighbor_info.append({
                    'guapos_index': int(guapos_idx),
                    'neighbor_index': int(pubchem_index),
                    'distance': float(distances[i][j]),
                    'guapos_point': [float(guapos_point[0]), float(guapos_point[1])],
                    'neighbor_point': [float(pubchem_points[neighbor_pos][0]), float(pubchem_points[neighbor_pos][1])],
                    'guapos_label': guapos_label,
                    'neighbor_label': neighbor_label,
                    'guapos_formula': guapos_formula,
                    'neighbor_formula': neighbor_formula
                })
    
    return neighbor_info

def create_interactive_plot(plot_data, knn_neighbors=None):
    """Create interactive Plotly visualization"""
    
    # Create DataFrame from plot data
    df_points = pd.DataFrame({
        'x': [point[0] for point in plot_data['points']],
        'y': [point[1] for point in plot_data['points']],
        'label': plot_data['labels'],
        'database': plot_data['databases'],
        'formula': plot_data['formulas'],
        'color': plot_data['colors']
    })
    
    # Create the main scatter plot
    fig = go.Figure()
    
    # Add points for each database
    databases = df_points['database'].unique()
    for db in databases:
        db_mask = df_points['database'] == db
        db_data = df_points[db_mask]
        
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
    
    # Add KNN neighbors from PubChem (red points) - with validation
    if knn_neighbors:
        # Filtrar índices válidos
        valid_neighbors = []
        for neighbor in knn_neighbors:
            guapos_idx = neighbor.get('guapos_index')
            neighbor_idx = neighbor.get('neighbor_index')
            
            # Verificar que ambos índices sean válidos
            if (guapos_idx is not None and neighbor_idx is not None and
                guapos_idx < len(plot_data['points']) and 
                neighbor_idx < len(plot_data['points'])):
                valid_neighbors.append(neighbor)
        
        if valid_neighbors:
            neighbor_indices = [n['neighbor_index'] for n in valid_neighbors]
            neighbor_points = df_points.iloc[neighbor_indices]
            
            # Create hover text for neighbors
            neighbor_hover_text = []
            for _, row in neighbor_points.iterrows():
                text = f"Database: KNN Neighbor (PubChem)<br>"
                if row['label']:
                    text += f"Name: {row['label']}<br>"
                if row['formula']:
                    formatted_formula = format_chemical_formula(row['formula'])
                    text += f"Formula: {formatted_formula}"
                neighbor_hover_text.append(text)
            
            fig.add_trace(go.Scatter(
                x=neighbor_points['x'],
                y=neighbor_points['y'],
                mode='markers',
                marker=dict(
                    color='red',
                    size=10,
                    opacity=0.7
                ),
                name='KNN Neighbors',
                hovertext=neighbor_hover_text,
                hoverinfo='text',
                showlegend=True
            ))
            
            # Add connections between GUAPOS points and their neighbors
            for neighbor in valid_neighbors:
                guapos_idx = neighbor['guapos_index']
                neighbor_idx = neighbor['neighbor_index']
                
                x0, y0 = plot_data['points'][guapos_idx]
                x1, y1 = plot_data['points'][neighbor_idx]
                
                fig.add_trace(go.Scatter(
                    x=[x0, x1],
                    y=[y0, y1],
                    mode='lines',
                    line=dict(color='red', width=1, dash='solid'),
                    opacity=0.6,
                    showlegend=False,
                    hoverinfo='skip'
                ))
    
    # Add KNN connections if they exist (from precomputed data) - with validation
    if plot_data.get('knn_connections'):
        connections = plot_data['knn_connections']
        
        # Filtrar conexiones válidas
        valid_connections = []
        for conn in connections:
            guapos_idx = conn.get('guapos_index')
            neighbor_idx = conn.get('neighbor_index')
            
            # Verificar que ambos índices sean válidos
            if (guapos_idx is not None and neighbor_idx is not None and
                guapos_idx < len(plot_data['points']) and 
                neighbor_idx < len(plot_data['points'])):
                valid_connections.append(conn)
        
        if valid_connections:
            neighbor_indices = [conn['neighbor_index'] for conn in valid_connections]
            neighbor_points = df_points.iloc[neighbor_indices]
            
            # Create hover text for neighbors
            neighbor_hover_text = []
            for _, row in neighbor_points.iterrows():
                text = f"Database: KNN Neighbor<br>"
                if row['label']:
                    text += f"Name: {row['label']}<br>"
                if row['formula']:
                    formatted_formula = format_chemical_formula(row['formula'])
                    text += f"Formula: {formatted_formula}"
                neighbor_hover_text.append(text)
            
            # Add red points for KNN neighbors
            fig.add_trace(go.Scatter(
                x=neighbor_points['x'],
                y=neighbor_points['y'],
                mode='markers',
                marker=dict(
                    color='red',
                    size=10,
                    opacity=0.7
                ),
                name='KNN Neighbors',
                hovertext=neighbor_hover_text,
                hoverinfo='text',
                showlegend=True
            ))
            
            # Add connections between GUAPOS points and their neighbors
            for connection in valid_connections:
                guapos_idx = connection['guapos_index']
                neighbor_idx = connection['neighbor_index']
                
                x0, y0 = plot_data['points'][guapos_idx]
                x1, y1 = plot_data['points'][neighbor_idx]
                
                fig.add_trace(go.Scatter(
                    x=[x0, x1],
                    y=[y0, y1],
                    mode='lines',
                    line=dict(color='red', width=1, dash='solid'),
                    opacity=0.6,
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
        height=700,
        width=1000
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

def create_molecule_info_panel(plot_data, selected_point):
    """Create molecule information panel"""
    if selected_point is None:
        return None
    
    point_idx = selected_point['pointIndex']
    
    if point_idx >= len(plot_data['points']):
        return None
    
    info = {
        'Database': plot_data['databases'][point_idx],
        'Name': plot_data['labels'][point_idx] or 'Unknown',
        'Formula': plot_data['formulas'][point_idx] or 'Unknown',
        'Coordinates': f"({plot_data['points'][point_idx][0]:.3f}, {plot_data['points'][point_idx][1]:.3f})"
    }
    
    return info

def main():
    # Add the header image and title
    st.image("NGC6523_BVO_2.jpg", use_column_width=True)
    
    col1, col2 = st.columns([1, 3])
    with col1:
        st.empty()
        
    with col2:
        st.markdown('<p class="main-title">AI-ITACA | Artificial Intelligence Integral Tool for AstroChemical Analysis</p>', unsafe_allow_html=True)
    
    st.markdown("""
    A remarkable upsurge in the complexity of molecules identified in the interstellar medium (ISM) is currently occurring, with over 80 new species discovered in the last three years. A number of them have been emphasized by prebiotic experiments as vital molecular building blocks of life. Since our Solar System was formed from a molecular cloud in the ISM, it prompts the query as to whether the rich interstellar chemical reservoir could have played a role in the emergence of life. The improved sensitivities of state-of-the-art astronomical facilities, such as the Atacama Large Millimeter/submillimeter Array (ALMA) and the James Webb Space Telescope (JWST), are revolutionizing the discovery of new molecules in space. However, we are still just scraping the tip of the iceberg. We are far from knowing the complete catalogue of molecules that astrochemistry can offer, as well as the complexity they can reach.<br><br>
    <strong>Artificial Intelligence Integral Tool for AstroChemical Analysis (AI-ITACA)</strong>, proposes to combine complementary machine learning (ML) techniques to address all the challenges that astrochemistry is currently facing. AI-ITACA will significantly contribute to the development of new AI-based cutting-edge analysis software that will allow us to make a crucial leap in the characterization of the level of chemical complexity in the ISM, and in our understanding of the contribution that interstellar chemistry might have in the origin of life.
    """, unsafe_allow_html=True)
    
    st.markdown("""
    <div class="info-box">
    <h4>About GUAPOS</h4>
    <p>The G31.41+0.31 Unbiased ALMA sPectral Observational Survey (GUAPOS) project targets the hot molecular core (HMC) G31.41+0.31 (G31) to reveal the complex chemistry of one of the most chemically rich high-mass star-forming regions outside the Galactic center (GC).</p>
    </div>
    """, unsafe_allow_html=True)
    
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
        
        if uploaded_file is not None:
            try:
                plot_data = json.load(uploaded_file)
                st.success("✅ Data loaded successfully!")
            except Exception as e:
                st.error(f"❌ Error loading file: {str(e)}")
                plot_data = None
        else:
            # Try to load from default path
            try:
                plot_data = load_plot_data("plot_data.json")
                if plot_data:
                    st.info("📁 Using default plot_data.json file")
            except:
                plot_data = None
                st.warning("⚠️ Please upload a plot data JSON file")
        
        if plot_data:
            st.markdown("---")
            st.header("🎨 Display Options")
            
            # Show database summary
            db_counts = pd.Series(plot_data['databases']).value_counts()
            st.write("**Database Summary:**")
            for db, count in db_counts.items():
                st.write(f"• {db}: {count} molecules")
            
            # Show KNN info if available
            if plot_data.get('knn_connections'):
                st.write(f"• KNN Connections: {len(plot_data['knn_connections'])}")
            
            st.markdown("---")
            st.header("🔍 Interaction Guide")
            st.info("""
            - **Hover** over points to see molecule details
            - **Click** on points to view detailed information
            - **Zoom** with mouse wheel or touchpad
            - **Pan** by dragging the plot
            - **Reset view** with home button in toolbar
            """)
    
    # Main content
    if plot_data is None:
        st.warning("""
        ## Welcome to Chemical Space Visualization!
        
        To get started:
        
        1. Run the main analysis script to generate `plot_data.json`
        2. Upload the JSON file using the sidebar
        3. Or place `plot_data.json` in the same directory as this app
        
        The visualization will show:
        - Molecular distributions across different databases
        - PCA/UMAP projections of chemical space
        - KNN connections between molecules
        - Interactive exploration capabilities
        """)
        
        # Show example image or placeholder
        st.image("https://via.placeholder.com/800x400/1f77b4/ffffff?text=Chemical+Space+Visualization", 
                use_column_width=True)
        
        return
    
    # Find KNN neighbors
    knn_neighbors = find_knn_neighbors(plot_data)
    
    # Create tabs for different views
    tab1, tab2, tab3, tab4 = st.tabs(["📈 Interactive Map", "📊 Statistics", "🔍 Molecule Explorer", "🧲 KNN Neighbors"])
    
    with tab1:
        st.markdown('<h2 class="sub-header">Interactive Chemical Space Map</h2>', unsafe_allow_html=True)
        
        # Create interactive plot
        fig = create_interactive_plot(plot_data, knn_neighbors)
        
        # Display the plot
        st.plotly_chart(fig, use_container_width=True)
        
        # Add some context
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("""
            <div class="info-box">
            <h4>About this Visualization</h4>
            <p>This map shows molecules projected into 2D space using PCA and UMAP dimensionality reduction techniques. 
            Colors represent different databases, and connections show molecular similarities detected by KNN analysis.</p>
            </div>
            """, unsafe_allow_html=True)
        
        with col2:
            st.markdown("""
            <div class="info-box">
            <h4>How to Interpret</h4>
            <p>• <strong>Closer points</strong> = More similar molecules<br>
            • <strong>Colors</strong> = Source database<br>
            • <strong>Gray lines</strong> = KNN similarity connections<br>
            • <strong>Red points/lines</strong> = KNN neighbors from PubChem close to GUAPOS molecules<br>
            • <strong>Larger points</strong> = GUAPOS molecules</p>
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
            # Additional statistics
            total_molecules = len(plot_data['points'])
            unique_databases = len(set(plot_data['databases']))
            
            st.metric("Total Molecules", f"{total_molecules:,}")
            st.metric("Unique Databases", unique_databases)
            
            if plot_data.get('knn_connections'):
                st.metric("KNN Connections", len(plot_data['knn_connections']))
            
            if knn_neighbors:
                st.metric("KNN Neighbors Found", len(knn_neighbors))
            
            # Show database details
            st.write("**Database Details:**")
            db_counts = pd.Series(plot_data['databases']).value_counts()
            for db, count in db_counts.items():
                percentage = (count / total_molecules) * 100
                st.write(f"• {db}: {count} ({percentage:.1f}%)")
    
    with tab3:
        st.markdown('<h2 class="sub-header">Molecule Explorer</h2>', unsafe_allow_html=True)
        
        # Create a searchable table of molecules
        molecules_df = pd.DataFrame({
            'Database': plot_data['databases'],
            'Name': plot_data['labels'],
            'Formula': plot_data['formulas'],
            'X': [p[0] for p in plot_data['points']],
            'Y': [p[1] for p in plot_data['points']]
        })
        
        # Search and filter options
        col1, col2 = st.columns(2)
        
        with col1:
            selected_db = st.multiselect(
                "Filter by Database",
                options=sorted(molecules_df['Database'].unique()),
                default=sorted(molecules_df['Database'].unique())
            )
        
        with col2:
            search_term = st.text_input("Search by Name or Formula", "")
        
        # Apply filters
        filtered_df = molecules_df[molecules_df['Database'].isin(selected_db)]
        if search_term:
            filtered_df = filtered_df[
                filtered_df['Name'].str.contains(search_term, case=False, na=False) |
                filtered_df['Formula'].str.contains(search_term, case=False, na=False)
            ]
        
        # Display results
        st.write(f"**Found {len(filtered_df)} molecules**")
        
        # Pagination
        page_size = 20
        total_pages = max(1, (len(filtered_df) + page_size - 1) // page_size)
        page = st.number_input("Page", min_value=1, max_value=total_pages, value=1)
        
        start_idx = (page - 1) * page_size
        end_idx = min(start_idx + page_size, len(filtered_df))
        
        # Display molecules
        for idx in range(start_idx, end_idx):
            molecule = filtered_df.iloc[idx]
            with st.expander(f"{molecule['Name'] or 'Unknown'} - {molecule['Database']}"):
                col1, col2 = st.columns(2)
                
                with col1:
                    st.write(f"**Database:** {molecule['Database']}")
                    st.write(f"**Name:** {molecule['Name'] or 'Unknown'}")
                    st.write(f"**Formula:** {format_chemical_formula(molecule['Formula'])}", unsafe_allow_html=True)
                
                with col2:
                    st.write(f"**Coordinates:** ({molecule['X']:.3f}, {molecule['Y']:.3f})")
                    st.write(f"**Index:** {filtered_df.index[idx]}")
        
        # Download option
        csv = filtered_df.to_csv(index=False)
        st.download_button(
            label="📥 Download Filtered Data",
            data=csv,
            file_name="filtered_molecules.csv",
            mime="text/csv"
        )
    
    with tab4:
        st.markdown('<h2 class="sub-header">KNN Neighbors Analysis</h2>', unsafe_allow_html=True)
        
        # Use precomputed KNN connections if available, otherwise use the ones we found
        if plot_data.get('knn_connections'):
            knn_data = plot_data['knn_connections']
            st.info("Using precomputed KNN connections from the uploaded data.")
        else:
            knn_data = knn_neighbors
            st.info("Using KNN neighbors calculated in real-time.")
        
        if not knn_data:
            st.info("No KNN neighbors found. Make sure you have GUAPOS points and PubChem points in your data.")
        else:
            st.metric("Total KNN Neighbors Found", len(knn_data))
            
            # Create a summary table with validation
            neighbor_summary = []
            for neighbor in knn_data:
                # Handle both formats (precomputed and real-time)
                guapos_idx = neighbor.get('guapos_index')
                neighbor_idx = neighbor.get('neighbor_index')
                
                # Validar índices antes de acceder
                if (guapos_idx is not None and neighbor_idx is not None and
                    guapos_idx < len(plot_data['labels']) and 
                    neighbor_idx < len(plot_data['labels'])):
                    
                    guapos_name = plot_data['labels'][guapos_idx] or f"Index {guapos_idx}"
                    neighbor_name = plot_data['labels'][neighbor_idx] or f"Index {neighbor_idx}"
                    guapos_formula = plot_data['formulas'][guapos_idx] if guapos_idx < len(plot_data['formulas']) else 'Unknown'
                    neighbor_formula = plot_data['formulas'][neighbor_idx] if neighbor_idx < len(plot_data['formulas']) else 'Unknown'
                    distance = neighbor.get('distance', 'N/A')
                    
                    neighbor_summary.append({
                        'GUAPOS Molecule': guapos_name,
                        'GUAPOS Formula': guapos_formula,
                        'Neighbor Molecule': neighbor_name,
                        'Neighbor Formula': neighbor_formula,
                        'Distance': f"{distance:.4f}" if isinstance(distance, (int, float)) else distance
                    })
            
            if neighbor_summary:
                neighbor_df = pd.DataFrame(neighbor_summary)
                
                # Group by GUAPOS molecule to show count
                guapos_counts = neighbor_df['GUAPOS Molecule'].value_counts().reset_index()
                guapos_counts.columns = ['GUAPOS Molecule', 'Number of Neighbors']
                
                col1, col2 = st.columns(2)
                
                with col1:
                    st.write("**Neighbors per GUAPOS Molecule:**")
                    st.dataframe(guapos_counts, use_container_width=True)
                
                with col2:
                    st.write("**All KNN Neighbors:**")
                    st.dataframe(neighbor_df, use_container_width=True)
                
                # Download option
                csv = neighbor_df.to_csv(index=False)
                st.download_button(
                    label="📥 Download KNN Neighbors Data",
                    data=csv,
                    file_name="knn_neighbors.csv",
                    mime="text/csv"
                )
            else:
                st.warning("No valid KNN neighbors found after validation.")

if __name__ == "__main__":
    main()
