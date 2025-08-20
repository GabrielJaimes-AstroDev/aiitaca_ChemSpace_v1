import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import json
import re

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
        'index': range(len(plot_data['points']))
    })
    
    # Identify neighbor indices
    neighbor_indices = set()
    if plot_data.get('knn_connections'):
        for connection in plot_data['knn_connections']:
            neighbor_indices.add(connection['neighbor_index'])
    
    # Create the main scatter plot
    fig = go.Figure()
    
    # Add points for each database (excluding neighbors which will be highlighted separately)
    databases = df_points['database'].unique()
    for db in databases:
        db_mask = (df_points['database'] == db) & (~df_points['index'].isin(neighbor_indices))
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

def create_neighbor_analysis_panel(plot_data):
    """Create panel showing neighbor candidate analysis"""
    if not plot_data.get('knn_connections'):
        return None
    
    # Get unique neighbor indices
    neighbor_indices = set()
    for connection in plot_data['knn_connections']:
        neighbor_indices.add(connection['neighbor_index'])
    
    # Count neighbors by database
    neighbor_dbs = {}
    for idx in neighbor_indices:
        if idx < len(plot_data['databases']):
            db = plot_data['databases'][idx]
            neighbor_dbs[db] = neighbor_dbs.get(db, 0) + 1
    
    return neighbor_dbs

def create_molecule_info_panel(plot_data, selected_point):
    """Create molecule information panel"""
    if selected_point is None:
        return None
    
    point_idx = selected_point['pointIndex']
    
    if point_idx >= len(plot_data['points']):
        return None
    
    # Check if this is a neighbor candidate
    is_neighbor = False
    if plot_data.get('knn_connections'):
        neighbor_indices = set()
        for connection in plot_data['knn_connections']:
            neighbor_indices.add(connection['neighbor_index'])
        is_neighbor = point_idx in neighbor_indices
    
    info = {
        'Database': plot_data['databases'][point_idx],
        'Name': plot_data['labels'][point_idx] or 'Unknown',
        'Formula': plot_data['formulas'][point_idx] or 'Unknown',
        'Coordinates': f"({plot_data['points'][point_idx][0]:.3f}, {plot_data['points'][point_idx][1]:.3f})",
        'Is Neighbor Candidate': is_neighbor
    }
    
    return info

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
            
            # Show neighbor analysis if available
            neighbor_dbs = create_neighbor_analysis_panel(plot_data)
            if neighbor_dbs:
                st.write("**Neighbor Candidates by Database:**")
                for db, count in neighbor_dbs.items():
                    st.write(f"• {db}: {count} candidates")
            
            # Show KNN info if available
            if plot_data.get('knn_connections'):
                st.write(f"• KNN Connections: {len(plot_data['knn_connections'])}")
            
            st.markdown("---")
            st.header("🔍 Interaction Guide")
            st.info("""
            - **Hover** over points to see molecule details
            - **Red points** = Neighbor candidates
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
        - Neighbor candidates highlighted in red
        - Interactive exploration capabilities
        """)
        
        # Show example image or placeholder
        st.image("https://via.placeholder.com/800x400/1f77b4/ffffff?text=Chemical+Space+Visualization", 
                use_column_width=True)
        
        return
    
    # Create tabs for different views
    tab1, tab2, tab3, tab4 = st.tabs(["📈 Interactive Map", "📊 Statistics", "🔍 Molecule Explorer", "🎯 Neighbor Analysis"])
    
    with tab1:
        st.markdown('<h2 class="sub-header">Interactive Chemical Space Map</h2>', unsafe_allow_html=True)
        
        # Create interactive plot
        fig = create_interactive_plot(plot_data)
        
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
            • <strong>Red points</strong> = Neighbor candidates<br>
            • <strong>Gray lines</strong> = KNN similarity connections<br>
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
            
            # Show neighbor candidate count
            neighbor_dbs = create_neighbor_analysis_panel(plot_data)
            if neighbor_dbs:
                total_neighbors = sum(neighbor_dbs.values())
                st.metric("Neighbor Candidates", total_neighbors)
            
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
            'Y': [p[1] for p in plot_data['points']],
            'Index': range(len(plot_data['points']))
        })
        
        # Identify neighbor candidates
        neighbor_indices = set()
        if plot_data.get('knn_connections'):
            for connection in plot_data['knn_connections']:
                neighbor_indices.add(connection['neighbor_index'])
        
        molecules_df['Is_Neighbor'] = molecules_df['Index'].isin(neighbor_indices)
        
        # Search and filter options
        col1, col2, col3 = st.columns(3)
        
        with col1:
            selected_db = st.multiselect(
                "Filter by Database",
                options=sorted(molecules_df['Database'].unique()),
                default=sorted(molecules_df['Database'].unique())
            )
        
        with col2:
            neighbor_filter = st.selectbox(
                "Neighbor Status",
                options=["All", "Neighbors Only", "Non-Neighbors Only"]
            )
        
        with col3:
            search_term = st.text_input("Search by Name or Formula", "")
        
        # Apply filters
        filtered_df = molecules_df[molecules_df['Database'].isin(selected_db)]
        
        if neighbor_filter == "Neighbors Only":
            filtered_df = filtered_df[filtered_df['Is_Neighbor'] == True]
        elif neighbor_filter == "Non-Neighbors Only":
            filtered_df = filtered_df[filtered_df['Is_Neighbor'] == False]
        
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
            card_class = "neighbor-point" if molecule['Is_Neighbor'] else "molecule-card"
            
            with st.expander(f"{molecule['Name'] or 'Unknown'} - {molecule['Database']} {'🎯' if molecule['Is_Neighbor'] else ''}"):
                st.markdown(f"<div class='{card_class}'>", unsafe_allow_html=True)
                
                col1, col2 = st.columns(2)
                
                with col1:
                    st.write(f"**Database:** {molecule['Database']}")
                    st.write(f"**Name:** {molecule['Name'] or 'Unknown'}")
                    if molecule['Is_Neighbor']:
                        st.write("**Status:** 🎯 Neighbor Candidate")
                    
                with col2:
                    st.write(f"**Formula:** {format_chemical_formula(molecule['Formula'])}", unsafe_allow_html=True)
                    st.write(f"**Coordinates:** ({molecule['X']:.3f}, {molecule['Y']:.3f})")
                    st.write(f"**Index:** {molecule['Index']}")
                
                st.markdown("</div>", unsafe_allow_html=True)
        
        # Download option
        csv = filtered_df.to_csv(index=False)
        st.download_button(
            label="📥 Download Filtered Data",
            data=csv,
            file_name="filtered_molecules.csv",
            mime="text/csv"
        )
    
    with tab4:
        st.markdown('<h2 class="sub-header">🎯 Neighbor Candidate Analysis</h2>', unsafe_allow_html=True)
        
        if not plot_data.get('knn_connections'):
            st.info("No KNN connections found in the data.")
        else:
            # Get neighbor analysis
            neighbor_dbs = create_neighbor_analysis_panel(plot_data)
            
            if neighbor_dbs:
                col1, col2 = st.columns(2)
                
                with col1:
                    # Neighbor distribution by database
                    neighbor_df = pd.DataFrame({
                        'Database': list(neighbor_dbs.keys()),
                        'Count': list(neighbor_dbs.values())
                    })
                    
                    fig = px.pie(
                        neighbor_df,
                        values='Count',
                        names='Database',
                        title='Neighbor Candidates by Database',
                        color='Database',
                        color_discrete_map={
                            'GUAPOS': '#1f77b4',
                            'TMC_1': '#ff7f0e',
                            'All_Discoveries': '#2ca02c',
                            'KIDA': '#FFFF00',
                            'PubChem': '#9467bd'
                        }
                    )
                    st.plotly_chart(fig, use_container_width=True)
                
                with col2:
                    # Neighbor statistics
                    total_neighbors = sum(neighbor_dbs.values())
                    total_connections = len(plot_data['knn_connections'])
                    
                    st.metric("Total Neighbor Candidates", total_neighbors)
                    st.metric("Total KNN Connections", total_connections)
                    st.metric("Average Connections per Candidate", f"{total_connections/total_neighbors:.1f}")
                    
                    st.write("**Neighbor Distribution:**")
                    for db, count in neighbor_dbs.items():
                        percentage = (count / total_neighbors) * 100
                        st.write(f"• {db}: {count} ({percentage:.1f}%)")
                
                # Show detailed neighbor information
                st.markdown("### Detailed Neighbor Information")
                
                # Get all neighbor molecules
                neighbor_indices = set()
                for connection in plot_data['knn_connections']:
                    neighbor_indices.add(connection['neighbor_index'])
                
                neighbor_molecules = []
                for idx in neighbor_indices:
                    if idx < len(plot_data['points']):
                        neighbor_molecules.append({
                            'Index': idx,
                            'Database': plot_data['databases'][idx],
                            'Name': plot_data['labels'][idx] or 'Unknown',
                            'Formula': plot_data['formulas'][idx] or 'Unknown',
                            'X': plot_data['points'][idx][0],
                            'Y': plot_data['points'][idx][1]
                        })
                
                neighbor_df = pd.DataFrame(neighbor_molecules)
                
                if not neighbor_df.empty:
                    # Display neighbor table
                    st.dataframe(
                        neighbor_df,
                        column_config={
                            "Index": st.column_config.NumberColumn("Index"),
                            "Database": "Database",
                            "Name": "Name",
                            "Formula": "Formula",
                            "X": st.column_config.NumberColumn("X", format="%.3f"),
                            "Y": st.column_config.NumberColumn("Y", format="%.3f")
                        },
                        hide_index=True,
                        use_container_width=True
                    )
                    
                    # Download neighbors
                    csv = neighbor_df.to_csv(index=False)
                    st.download_button(
                        label="📥 Download Neighbor Candidates",
                        data=csv,
                        file_name="neighbor_candidates.csv",
                        mime="text/csv"
                    )
            else:
                st.info("No neighbor candidates found in the data.")

if __name__ == "__main__":
    main()
