<style>
.loader {
  border: 5px solid #f3f3f3; /* Light grey */
  border-top: 5px solid #3498db; /* Blue */
  border-radius: 50%;
  width: 50px;
  height: 50px;
  animation: spin 2s linear infinite;
  display: inline-block;
}

@keyframes spin {
  0% {
    transform: rotate(0deg);
  }
  100% {
    transform: rotate(360deg);
  }
}
.center {
  position:fixed;
  left: 50%;
  top:50%;
}
</style>
<template>
  <div>
    <div id="aligmentview" style="width:100%"></div>
  </div>
</template>
<script>
/* eslint-disable */
import { AlignmentViewer } from "@/amromicsjs";
import EventBus from "@/event-bus.js";
import SampleAPI from "@/api/SampleAPI";
// import SampleIGV from "@/components/Visualization/IGV";
export default {
  name: "AlignmentComp",
  props: ["alignmentData"],
  data() {
    return {
      loading: false,
      list_alignments: [],
      current_alignment:undefined,
      alignmentview:undefined
    };
  },
  computed: {
    collectionId() {
      return this.$route.params.cid;
      ;
    }
  },
  async mounted() {
    this.loading = true;
    var ctx = document.getElementById("aligmentview");
    //console.log(this.core_data);
    //console.log(Phylogeny);
    
    this.list_alignments = this.alignmentData.alignments;
    const value = await SampleAPI.fetchAlignment(this.collectionId,this.alignmentData.alignments[0].gene);
    this.current_alignment=value.data
    this.alignmentview = new AlignmentViewer(ctx);

    var tree_data = atob(this.alignmentData.alignments[0].tree).replace(
      /.ref/g,
      ""
    );
    tree_data = tree_data.replace(/.fasta/g, "");
    //console.log(this.alignmentData.alignments[0])
    this.alignmentview.load(
      this.alignmentData.alignments[0].gene,
      tree_data,
      this.current_alignment
    );
    //alignmentview.setOptions({width:ctx.clientWidth,height:0});
    
    this.alignmentview.draw();
   
    EventBus.$on("gene_id_emited", gene_id => {
      //console.log('gene_id_emited'+gene_id);
      this.loading=true;
      this.reloadAlignment(gene_id);
      
      this.loading=false;
     
    });
    EventBus.$on("samples_emited", arr_ids => {
      //console.log('sample_emited '+arr_ids);
      this.loading=true;
 
      this.alignmentview.setActiveNames(arr_ids);
      this.loading=false;
  
    });
    this.loading = false;
    //document.getElementById("al_loader").display="none";
  },
  async created() {},
  methods: {
    async reloadAlignment(gene_id){
      for (var i = 0; i < this.list_alignments.length; i++) {
        if (this.list_alignments[i].gene == gene_id) {
          var tree = atob(this.alignmentData.alignments[i].tree).replace(
            /.ref/g,
            ""
          );
          tree = tree.replace(/.fasta/g, "");
          //console.log(tree);
          const value = await SampleAPI.fetchAlignment(this.collectionId,gene_id);
          this.current_alignment=value.data
          this.alignmentview.load(
            this.alignmentData.alignments[i].gene,
            tree,
            this.current_alignment
          );
          this.alignmentview.draw();
          break;
        }
      }
    }
  }
};
</script>
