<style scoped>
</style>
<template>
<div id="browser" >
</div>
</template>
<script>
/* eslint-disable */
import {Browser} from "@/amromicsjs";
import EventBus from '@/event-bus.js';
// import SampleIGV from "@/components/Visualization/IGV";
export default {
  name: 'GenomeBrowser',
  props: ["list_contig","knowngene","GC_content","GC_skew"],
  data() {
    return {
      loading: false,
      current_zoom: 3,
      current_type: "prot",
      pos_track: 0,
      svg_title: undefined,
      scale:0.1,
      current_contig:undefined,
      selected_contig:undefined
    };
  },
  mounted() {
    this.loading = true;
    var ctx=document.getElementById('browser')
    //console.log(this.GC_content);
    var browser = new Browser(ctx);
    browser.load(this.list_contig,this.knowngene,this.GC_skew,this.GC_content);
    browser.setOptions({width:ctx.clientWidth});
    browser.draw();
    EventBus.$on('contig_emited', contig_id => {
      
      browser.changeContig(contig_id);
    });
    EventBus.$on('element_emited', element => {
      //console.log(element);
      browser.changeContig(element.contig);
      browser.locatePosition(element.pos);
    });
    this.loading = false;
  },
  async created() {


  },
  method: {

}


};
</script>
